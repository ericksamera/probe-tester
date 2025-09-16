"""Helpers for retrieving CSL JSON and BibTeX records for DOIs.

The functions in this module centralise the logic for making HTTP requests to
``doi.org`` and the Crossref API.  They normalise DOI strings, attach a
meaningful ``User-Agent`` header (required by Crossref) and provide a
consistent :class:`CitationFetchError` exception when a fetch fails.  This
allows calling code to gracefully handle transient network problems instead of
propagating raw ``URLError``/``HTTPError`` exceptions.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any, Callable, Iterable, List
from urllib import error, parse, request

__all__ = [
    "CitationFetchError",
    "DEFAULT_TIMEOUT",
    "fetch_csl_via_doi",
    "fetch_csl_via_crossref",
    "fetch_bibtex_via_doi",
    "fetch_bibtex_via_crossref",
    "fetch_csl",
    "fetch_bibtex",
    "normalize_doi",
]

DEFAULT_TIMEOUT = 10.0

# Crossref asks for a descriptive user-agent containing contact details.  The
# default value is intentionally generic, but the environment variable allows
# deployments to provide a more descriptive string without modifying code.
DEFAULT_USER_AGENT = (
    "PaperClip/1.0 (+https://example.com/; mailto:paperclip@example.com)"
)


def _user_agent() -> str:
    """Return the user-agent string to send with outbound HTTP requests."""

    return os.environ.get("PAPERCLIP_USER_AGENT", DEFAULT_USER_AGENT)


class CitationFetchError(RuntimeError):
    """Raised when a DOI reference cannot be retrieved."""


_DOI_PREFIX_RE = re.compile(r"^(?:urn:)?doi:\s*", re.IGNORECASE)
_DOI_LEADING_WRAPPERS = {
    "<": ">",
    "[": "]",
    "{": "}",
    "(": ")",
    "\"": "\"",
    "'": "'",
}
_DOI_HOSTS = {"doi.org", "dx.doi.org"}
_TRAILING_PUNCTUATION = ",.;"


def _strip_wrapping(value: str) -> str:
    """Remove matching leading/trailing wrapper characters from ``value``."""

    while value and value[0] in _DOI_LEADING_WRAPPERS:
        closing = _DOI_LEADING_WRAPPERS[value[0]]
        if value.endswith(closing):
            value = value[1:-1].strip()
        else:
            break
    return value


def normalize_doi(doi: str) -> str:
    """Return the DOI identifier component of ``doi``.

    The helper accepts common DOI representations including raw identifiers,
    ``doi:``/``urn:doi:`` prefixes, wrapped values (``<...>``) and URLs pointing
    at ``doi.org``/``dx.doi.org``.  Any query strings or fragments are removed
    because they are not part of the DOI itself.  The result is percent-decoded
    and stripped of incidental punctuation such as trailing commas.
    """

    if not isinstance(doi, str):
        raise CitationFetchError("DOI must be provided as a string")

    cleaned = doi.strip()
    if not cleaned:
        raise CitationFetchError("DOI must be a non-empty string")

    cleaned = _strip_wrapping(cleaned)
    cleaned = _DOI_PREFIX_RE.sub("", cleaned, count=1)

    lowered = cleaned.lower()
    if lowered.startswith("doi="):
        cleaned = cleaned[4:]
        lowered = cleaned.lower()

    if lowered.startswith("https://") or lowered.startswith("http://"):
        parsed = parse.urlparse(cleaned)
        host = parsed.netloc.lower()
        if host in _DOI_HOSTS or any(host.endswith(f".{known}") for known in _DOI_HOSTS):
            cleaned = parsed.path.lstrip("/")
        else:
            cleaned = parsed.geturl()
        if parsed.query:
            cleaned = cleaned.split("?", 1)[0]
        if parsed.fragment:
            cleaned = cleaned.split("#", 1)[0]
    elif lowered.startswith("doi.org/") or lowered.startswith("dx.doi.org/"):
        cleaned = cleaned.split("/", 1)[1]
        cleaned = cleaned.split("?", 1)[0]
        cleaned = cleaned.split("#", 1)[0]

    cleaned = cleaned.strip()
    if not cleaned:
        raise CitationFetchError("DOI must be a non-empty string")

    cleaned = parse.unquote(cleaned)
    cleaned = re.sub(r"\s+", "", cleaned)
    cleaned = cleaned.rstrip(_TRAILING_PUNCTUATION)
    cleaned = _strip_wrapping(cleaned)
    cleaned = cleaned.strip()

    if not cleaned:
        raise CitationFetchError("DOI must be a non-empty string")

    # Guard against stray control characters that can appear when copying from
    # PDFs or poorly formatted sources.
    if any(ord(ch) < 32 for ch in cleaned):
        raise CitationFetchError("DOI contains invalid control characters")

    return cleaned


def _build_doi_url(normalized_doi: str) -> str:
    return f"https://doi.org/{parse.quote(normalized_doi, safe='/')}"


def _build_crossref_url(normalized_doi: str, content_type: str) -> str:
    safe_doi = parse.quote(normalized_doi, safe="/")
    return f"https://api.crossref.org/works/{safe_doi}/transform/{content_type}"


@dataclass
class _FetchAttempt:
    description: str
    url: str
    accept: str
    parser: Callable[[str], Any]


def _http_get(url: str, accept: str, timeout: float) -> str:
    headers = {
        "Accept": accept,
        "User-Agent": _user_agent(),
    }
    req = request.Request(url, headers=headers)
    try:
        with request.urlopen(req, timeout=timeout) as response:
            charset = response.headers.get_content_charset() or "utf-8"
            data = response.read()
            try:
                return data.decode(charset)
            except UnicodeDecodeError as exc:  # pragma: no cover - validated in tests via mocks
                raise CitationFetchError(
                    f"Unable to decode response from {url} using charset '{charset}'"
                ) from exc
    except error.HTTPError as exc:  # pragma: no cover - exercised in tests via mocks
        raise CitationFetchError(
            f"HTTP error {exc.code} while fetching {url}: {exc.reason}"
        ) from exc
    except error.URLError as exc:  # pragma: no cover - exercised in tests via mocks
        raise CitationFetchError(f"Failed to fetch {url}: {exc.reason}") from exc
    except ValueError as exc:  # pragma: no cover - exercised in tests via mocks
        raise CitationFetchError(f"Invalid URL '{url}': {exc}") from exc


def _fetch_with_attempts(
    attempts: Iterable[_FetchAttempt],
    timeout: float,
    doi_for_error: str,
) -> Any:
    errors: List[str] = []
    for attempt in attempts:
        try:
            payload = _http_get(attempt.url, attempt.accept, timeout)
            return attempt.parser(payload)
        except CitationFetchError as exc:
            errors.append(f"{attempt.description}: {exc}")
    joined = "; ".join(errors)
    raise CitationFetchError(
        f"Unable to fetch data for DOI '{doi_for_error}': {joined}"
    )


def fetch_csl_via_doi(doi: str, timeout: float = DEFAULT_TIMEOUT) -> dict[str, Any]:
    """Fetch CSL JSON from ``doi.org`` for ``doi``."""

    normalized = normalize_doi(doi)

    def parse_json(payload: str) -> dict[str, Any]:
        try:
            return json.loads(payload)
        except json.JSONDecodeError as exc:  # pragma: no cover - validated via mocks
            raise CitationFetchError("Received invalid JSON from doi.org") from exc

    attempt = _FetchAttempt(
        description="doi.org CSL",
        url=_build_doi_url(normalized),
        accept="application/vnd.citationstyles.csl+json",
        parser=parse_json,
    )
    return _fetch_with_attempts([attempt], timeout, normalized)


def fetch_csl_via_crossref(doi: str, timeout: float = DEFAULT_TIMEOUT) -> dict[str, Any]:
    """Fetch CSL JSON by calling the Crossref ``transform`` endpoint."""

    normalized = normalize_doi(doi)

    def parse_json(payload: str) -> dict[str, Any]:
        try:
            return json.loads(payload)
        except json.JSONDecodeError as exc:  # pragma: no cover - validated via mocks
            raise CitationFetchError("Received invalid JSON from Crossref") from exc

    attempt = _FetchAttempt(
        description="Crossref CSL",
        url=_build_crossref_url(normalized, "application/vnd.citationstyles.csl+json"),
        accept="application/vnd.citationstyles.csl+json",
        parser=parse_json,
    )
    return _fetch_with_attempts([attempt], timeout, normalized)


def fetch_bibtex_via_doi(doi: str, timeout: float = DEFAULT_TIMEOUT) -> str:
    """Fetch a BibTeX entry from ``doi.org`` for ``doi``."""

    normalized = normalize_doi(doi)
    attempt = _FetchAttempt(
        description="doi.org BibTeX",
        url=_build_doi_url(normalized),
        accept="application/x-bibtex",
        parser=lambda payload: payload.strip(),
    )
    return _fetch_with_attempts([attempt], timeout, normalized)


def fetch_bibtex_via_crossref(doi: str, timeout: float = DEFAULT_TIMEOUT) -> str:
    """Fetch BibTeX from Crossref's ``transform`` endpoint."""

    normalized = normalize_doi(doi)
    attempt = _FetchAttempt(
        description="Crossref BibTeX",
        url=_build_crossref_url(normalized, "application/x-bibtex"),
        accept="application/x-bibtex",
        parser=lambda payload: payload.strip(),
    )
    return _fetch_with_attempts([attempt], timeout, normalized)


def fetch_csl(doi: str, timeout: float = DEFAULT_TIMEOUT) -> dict[str, Any]:
    """Fetch CSL JSON using doi.org first with a Crossref fallback."""

    normalized = normalize_doi(doi)
    attempts = [
        _FetchAttempt(
            description="doi.org CSL",
            url=_build_doi_url(normalized),
            accept="application/vnd.citationstyles.csl+json",
            parser=lambda payload: json.loads(payload),
        ),
        _FetchAttempt(
            description="Crossref CSL",
            url=_build_crossref_url(normalized, "application/vnd.citationstyles.csl+json"),
            accept="application/vnd.citationstyles.csl+json",
            parser=lambda payload: json.loads(payload),
        ),
    ]

    def parse_json(payload: str) -> dict[str, Any]:
        try:
            return json.loads(payload)
        except json.JSONDecodeError as exc:
            raise CitationFetchError("Received invalid JSON data") from exc

    # Replace parser to ensure consistent error handling while reusing URLs/accept
    attempts_with_parser = [
        _FetchAttempt(a.description, a.url, a.accept, parse_json) for a in attempts
    ]
    return _fetch_with_attempts(attempts_with_parser, timeout, normalized)


def fetch_bibtex(doi: str, timeout: float = DEFAULT_TIMEOUT) -> str:
    """Fetch BibTeX using doi.org first with a Crossref fallback."""

    normalized = normalize_doi(doi)
    attempts = [
        _FetchAttempt(
            description="doi.org BibTeX",
            url=_build_doi_url(normalized),
            accept="application/x-bibtex",
            parser=lambda payload: payload.strip(),
        ),
        _FetchAttempt(
            description="Crossref BibTeX",
            url=_build_crossref_url(normalized, "application/x-bibtex"),
            accept="application/x-bibtex",
            parser=lambda payload: payload.strip(),
        ),
    ]
    return _fetch_with_attempts(attempts, timeout, normalized)
