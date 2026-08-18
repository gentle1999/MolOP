"""Resolve a documentation deployment ref into mike version metadata."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass


_RELEASE_TAG = re.compile(
    r"^v(?P<version>[0-9]+\.[0-9]+\.[0-9]+"
    r"(?:(?:a|b|rc)[0-9]+|dev[0-9]*|\.(?:post|dev)[0-9]+)?)$"
)


@dataclass(frozen=True)
class DocsVersion:
    snapshot_kind: str
    version: str
    source_version: str
    channel: str
    update_latest: bool

    def github_output(self) -> str:
        values = {
            "snapshot_kind": self.snapshot_kind,
            "version": self.version,
            "source_version": self.source_version,
            "channel": self.channel,
            "update_latest": str(self.update_latest).lower(),
        }
        return "\n".join(f"{key}={value}" for key, value in values.items())


def resolve_docs_version(source_ref: str, *, update_latest: bool = False) -> DocsVersion:
    if source_ref == "main":
        return DocsVersion(
            snapshot_kind="development",
            version="main",
            source_version="",
            channel="development",
            update_latest=False,
        )

    match = _RELEASE_TAG.fullmatch(source_ref)
    if match is None:
        raise ValueError(
            "source_ref must be main or a PEP 440 release tag such as "
            f"v1.2.3, v1.2.3rc1, or v1.2.3.post1; got: {source_ref}"
        )

    version = match.group("version")
    return DocsVersion(
        snapshot_kind="tag",
        version=version,
        source_version=version,
        channel="release",
        update_latest=update_latest,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_ref")
    parser.add_argument("update_latest", choices=("true", "false"))
    args = parser.parse_args()
    print(
        resolve_docs_version(
            args.source_ref,
            update_latest=args.update_latest == "true",
        ).github_output()
    )


if __name__ == "__main__":
    main()
