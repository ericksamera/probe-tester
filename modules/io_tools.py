# modules/io_tools.py
import time
from typing import List
import subprocess
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def run_command(
        command: List[str],
        capture_output: bool = False,
        dry_run: bool = False
        ) -> subprocess.CompletedProcess[bytes] | None:
    """
    """
    logger.debug(f"Running: {' '.join(command)}")

    if dry_run:
        logger.debug("Dry-run mode: command not executed.")
        return None

    start = time.time()
    result = subprocess.run(command, capture_output=capture_output, check=True)
    end = time.time()

    duration = end - start
    logger.debug(f"Command completed in {duration:.2f}s")
    
    return result if capture_output else None

def ask(
        prompt: str,
        default: str | None = None
        ) -> str:
    """
    """

    show_default = f" [{default}]" if default else ""
    response = input(f"{prompt}{show_default}: ").strip()
    return response or (default if default is not None else "")

# ---