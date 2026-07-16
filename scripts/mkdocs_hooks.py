from __future__ import annotations

import warnings
from pathlib import Path

from mkdocs import plugins
from mkdocs_jupyter.plugin import NotebookFile
from nbformat.warnings import MissingIDFieldWarning


warnings.filterwarnings("ignore", category=MissingIDFieldWarning)

NOTEBOOK_NAMES = (
    "01-gaussian-parse-and-inspect",
    "02-batch-summary-filter-select",
    "03-transform-and-export",
)


@plugins.event_priority(-200)
def on_files(files, config):
    for file in list(files):
        normalized_src_uri = getattr(file, "norm_src_uri", file.src_uri)
        if normalized_src_uri.startswith("examples/") and normalized_src_uri.endswith(".ipynb"):
            files.remove(file)
            notebook_file = NotebookFile(
                file,
                use_directory_urls=config.use_directory_urls,
                site_dir=config.site_dir,
            )
            # i18n already removed the default-language source prefix from the
            # destination. Keep that mapping instead of deriving it from src_uri.
            notebook_file.dest_uri = file.dest_uri
            notebook_file.abs_dest_path = file.abs_dest_path
            notebook_file.url = file.url
            files.append(notebook_file)
    return files


@plugins.event_priority(-200)
def on_post_build(config):
    i18n = config.plugins.get("i18n")
    locale_prefixes = ("",)
    if i18n is not None:
        locale_prefixes = tuple(
            "" if locale == i18n.default_language else f"{locale}/"
            for locale in i18n.build_languages
        )

    for locale_prefix in locale_prefixes:
        examples_dir = Path(config.site_dir) / locale_prefix / "examples"
        for notebook_name in NOTEBOOK_NAMES:
            output_path = examples_dir / notebook_name / "index.html"
            if not output_path.is_file():
                raise RuntimeError(f"Notebook page was not built: {output_path}")
            output = output_path.read_text(encoding="utf-8")
            if not output.lstrip().lower().startswith("<!doctype html"):
                raise RuntimeError(f"Notebook page is not rendered HTML: {output_path}")
            if 'class="jp-Notebook"' not in output or 'class="jp-Cell' not in output:
                raise RuntimeError(f"Notebook cells are missing from page: {output_path}")
