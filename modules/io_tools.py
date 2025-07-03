"""
modules/io_tools.py

Utility functions for running shell commands and interacting with the user.
"""

import time
from typing import List, Optional, Any
import subprocess
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def run_command(
    command: List[str],
    capture_output: bool = False,
    dry_run: bool = False,
    check: bool = True,
    cwd: Optional[Path] = None,
    env: Optional[dict[str, str]] = None
) -> Optional[subprocess.CompletedProcess[Any]]:
    """
    Runs a shell command.

    Args:
        command: List of command-line arguments.
        capture_output: If True, capture stdout and stderr.
        dry_run: If True, only log the command, do not execute.
        check: If True, raises CalledProcessError on nonzero exit.
        cwd: Directory to execute the command in.
        env: Environment variables to use.

    Returns:
        CompletedProcess if capture_output is True, else None.

    Raises:
        subprocess.CalledProcessError if check is True and the command fails.
    """
    logger.debug("Running command: %s", " ".join(command))
    if dry_run:
        logger.info("Dry-run mode: command not executed.")
        return None

    start = time.time()
    try:
        result = subprocess.run(
            command,
            capture_output=capture_output,
            check=check,
            cwd=str(cwd) if cwd else None,
            env=env
        )
    except Exception as e:
        logger.error("Command failed: %s", e)
        raise
    end = time.time()

    logger.debug("Command completed in %.2fs", end - start)
    return result if capture_output else None

def ask(
    prompt: str,
    default: Optional[str] = None
) -> str:
    """
    Prompts the user for input, showing a default value.

    Args:
        prompt: The prompt text.
        default: Default value to use if user presses enter.

    Returns:
        User input, or the default value if input is empty.
    """
    show_default = f" [{default}]" if default else ""
    response = input(f"{prompt}{show_default}: ").strip()
    return response if response else (default if default is not None else "")

# ---