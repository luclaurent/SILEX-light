#!/usr/bin/env python3
"""Prepare Sphinx markdown output for publication on GitHub Wiki.

The Sphinx markdown builder produces files that are readable as plain markdown,
but GitHub Wiki benefits from a small amount of post-processing:

- publish only markdown pages, not builder internals such as ``.doctrees``
- rename ``index.md`` to ``Home.md`` for the wiki landing page
- rewrite internal ``.md`` links to wiki page links to avoid raw-file rendering
- preserve Sphinx API cross-links by injecting explicit anchors for function pages
"""

from __future__ import annotations

import argparse
import html
import posixpath
import re
import shutil
from pathlib import Path


MARKDOWN_LINK_RE = re.compile(r"(?P<prefix>\[[^\]]+\]\()(?P<target>[^)#]+\.md)(?P<suffix>#[^)]+)?\)")
API_HEADING_RE = re.compile(r"^(?P<level>#{3,6})\s+(?P<name>[A-Za-z0-9_.]+)\(.*\)\s*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Directory containing Sphinx markdown output")
    parser.add_argument("destination", type=Path, help="Directory for GitHub Wiki ready markdown")
    return parser.parse_args()


def destination_relative_path(source_relative_path: Path) -> Path:
    if source_relative_path == Path("index.md"):
        return Path("Home.md")
    return source_relative_path


def link_target_to_wiki_path(target: str, current_file: Path) -> str:
    target_path = Path(target)
    combined = Path(posixpath.normpath((current_file.parent / target_path).as_posix()))
    if combined.name == "index.md":
        combined = combined.with_name("Home.md")
    return combined.as_posix().removesuffix(".md")


def rewrite_links(markdown_text: str, current_relative_path: Path) -> str:
    def replacement(match: re.Match[str]) -> str:
        target = match.group("target")
        suffix = match.group("suffix") or ""
        wiki_target = link_target_to_wiki_path(target, current_relative_path)
        return f"{match.group('prefix')}{wiki_target}{suffix})"

    return MARKDOWN_LINK_RE.sub(replacement, markdown_text)


def inject_api_anchors(markdown_text: str) -> str:
    processed_lines: list[str] = []

    for line in markdown_text.splitlines():
        match = API_HEADING_RE.match(line)
        if match:
            anchor_name = html.escape(match.group("name"), quote=True)
            processed_lines.append(f'<a id="{anchor_name}"></a>')
        processed_lines.append(line)

    return "\n".join(processed_lines) + "\n"


def transform_markdown(markdown_text: str, current_relative_path: Path) -> str:
    transformed = rewrite_links(markdown_text, current_relative_path)
    transformed = inject_api_anchors(transformed)
    return transformed


def main() -> int:
    args = parse_args()
    source = args.source.resolve()
    destination = args.destination.resolve()

    if not source.is_dir():
        raise SystemExit(f"Source directory does not exist: {source}")

    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)

    markdown_files = sorted(source.rglob("*.md"))
    for source_file in markdown_files:
        source_relative_path = source_file.relative_to(source)
        destination_relative = destination_relative_path(source_relative_path)
        destination_file = destination / destination_relative
        destination_file.parent.mkdir(parents=True, exist_ok=True)

        source_text = source_file.read_text(encoding="utf-8")
        transformed_text = transform_markdown(source_text, source_relative_path)
        destination_file.write_text(transformed_text, encoding="utf-8")

    home_file = destination / "Home.md"
    if not home_file.exists():
        fallback_source = destination / "index.md"
        if fallback_source.exists():
            fallback_source.rename(home_file)
        else:
            home_file.write_text("# SILEXlight Wiki\n", encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())