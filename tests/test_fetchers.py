from __future__ import annotations

import json
import sys
from email.message import Message
from pathlib import Path
from typing import Iterator

import pytest

PROJECT_ROOT = str(Path(__file__).resolve().parents[1])
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from paperclip.fetchers import (
    CitationFetchError,
    _http_get,  # type: ignore[attr-defined]
    fetch_bibtex,
    fetch_csl,
    fetch_csl_via_doi,
    normalize_doi,
)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("10.1000/xyz", "10.1000/xyz"),
        (" DOI:10.1000/xyz ", "10.1000/xyz"),
        ("urn:doi:10.1000/xyz", "10.1000/xyz"),
        ("https://doi.org/10.1000/xyz", "10.1000/xyz"),
        ("https://dx.doi.org/10.1000/xyz", "10.1000/xyz"),
        ("doi.org/10.1000/xyz", "10.1000/xyz"),
        ("DX.DOI.ORG/10.1000/XYZ", "10.1000/XYZ"),
        ("<https://doi.org/10.1000/xyz>", "10.1000/xyz"),
        ("https://doi.org/10.1000/xyz?param=1", "10.1000/xyz"),
        ("https://doi.org/10.1000/xyz#fragment", "10.1000/xyz"),
        ("https://doi.org/10.1000/xyz%2Fabc", "10.1000/xyz/abc"),
        ("doi=10.1000/xyz.", "10.1000/xyz"),
        ("10.1000/xyz,", "10.1000/xyz"),
        ("\n10.1000/xyz\n", "10.1000/xyz"),
        ("<10.1000/xyz>;", "10.1000/xyz"),
        ("10.1000 / xyz", "10.1000/xyz"),
    ],
)
def test_normalize_doi_variants(value: str, expected: str) -> None:
    assert normalize_doi(value) == expected


@pytest.mark.parametrize("value", ["", "   ", "\t\n\r"])
def test_normalize_doi_rejects_empty(value: str) -> None:
    with pytest.raises(CitationFetchError):
        normalize_doi(value)


def test_normalize_doi_rejects_non_string() -> None:
    with pytest.raises(CitationFetchError):
        normalize_doi(123)  # type: ignore[arg-type]


def test_normalize_doi_rejects_control_characters() -> None:
    with pytest.raises(CitationFetchError):
        normalize_doi("10.1000/ab\x01c")


def _make_dummy_response(payload: bytes, content_type: str) -> object:
    class DummyResponse:
        def __init__(self, payload: bytes, content_type: str) -> None:
            self._payload = payload
            self.headers = Message()
            if content_type:
                self.headers.add_header("Content-Type", content_type)

        def read(self) -> bytes:
            return self._payload

        def __enter__(self) -> "DummyResponse":
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - nothing to clean up
            return None

    return DummyResponse(payload, content_type)


def test_fetch_uses_environment_user_agent(monkeypatch: pytest.MonkeyPatch) -> None:
    payload = json.dumps({"title": ["Example"]}).encode("utf-8")

    captured_headers: dict[str, str] = {}

    def fake_urlopen(req: object, timeout: float) -> object:
        # ``urllib.request.Request`` provides ``header_items`` for inspection.
        header_dict = {k.lower(): v for k, v in req.header_items()}  # type: ignore[attr-defined]
        captured_headers.update(header_dict)
        return _make_dummy_response(payload, "application/json; charset=utf-8")

    monkeypatch.setenv("PAPERCLIP_USER_AGENT", "CustomAgent/1.0")
    monkeypatch.setattr("paperclip.fetchers.request.urlopen", fake_urlopen)

    data = fetch_csl_via_doi("10.1000/xyz")

    assert data["title"] == ["Example"]
    assert captured_headers.get("user-agent") == "CustomAgent/1.0"


def test_http_get_unicode_decode_error(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_urlopen(req: object, timeout: float) -> object:
        return _make_dummy_response(b"\xff\xfe", "text/plain; charset=ascii")

    monkeypatch.setattr("paperclip.fetchers.request.urlopen", fake_urlopen)

    with pytest.raises(CitationFetchError) as excinfo:
        _http_get("https://example.invalid", "text/plain", timeout=1.0)

    assert "Unable to decode response" in str(excinfo.value)


def test_fetch_bibtex_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    def fake_http_get(url: str, accept: str, timeout: float) -> str:
        calls.append(url)
        if len(calls) == 1:
            raise CitationFetchError("network error")
        return "@article{foo, title={Example}}"

    monkeypatch.setattr("paperclip.fetchers._http_get", fake_http_get)

    result = fetch_bibtex("10.2000/example")

    assert result.startswith("@article")
    assert len(calls) == 2
    assert calls[0].startswith("https://doi.org/")
    assert calls[1].startswith("https://api.crossref.org/")


def test_fetch_bibtex_raises_combined_error(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_http_get(url: str, accept: str, timeout: float) -> str:
        raise CitationFetchError(f"boom from {url}")

    monkeypatch.setattr("paperclip.fetchers._http_get", fake_http_get)

    with pytest.raises(CitationFetchError) as excinfo:
        fetch_bibtex("10.9999/nowhere")

    message = str(excinfo.value)
    assert "Unable to fetch data for DOI '10.9999/nowhere'" in message
    assert "doi.org BibTeX" in message
    assert "Crossref BibTeX" in message


def test_fetch_csl_invalid_json_falls_back(monkeypatch: pytest.MonkeyPatch) -> None:
    responses: Iterator[str] = iter(["not json", json.dumps({"id": "ok"})])

    def fake_http_get(url: str, accept: str, timeout: float) -> str:
        return next(responses)

    monkeypatch.setattr("paperclip.fetchers._http_get", fake_http_get)

    data = fetch_csl("10.4000/example")

    assert data == {"id": "ok"}
