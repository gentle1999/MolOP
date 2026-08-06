from __future__ import annotations

import json
import re
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
NOTEBOOK_OUTPUT_RE = re.compile(
    r"^(?P<indent>[ \t]*)<!--\s*notebook-output:\s*"
    r"(?P<notebook>[^#\s]+)#(?P<cell_id>[A-Za-z0-9_.-]+)\s*-->[ \t]*$",
    re.MULTILINE,
)
PANDAS_STYLE_RE = re.compile(r"<style scoped>.*?</style>", re.DOTALL)


def _notebook_cell_html(notebook_path: Path, cell_id: str) -> str:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    cells = [cell for cell in notebook.get("cells", []) if cell.get("id") == cell_id]
    if len(cells) != 1:
        raise RuntimeError(
            f"Expected one notebook cell {cell_id!r} in {notebook_path}, found {len(cells)}"
        )

    fragments: list[str] = []
    for output in cells[0].get("outputs", []):
        html = output.get("data", {}).get("text/html")
        if isinstance(html, list):
            fragments.append("".join(html))
        elif isinstance(html, str):
            fragments.append(html)
    if not fragments:
        raise RuntimeError(
            f"Notebook cell {cell_id!r} in {notebook_path} has no saved text/html output. "
            "Execute the documentation notebooks before building MkDocs."
        )
    return PANDAS_STYLE_RE.sub("", "\n".join(fragments)).strip()


def on_page_markdown(markdown, page, config, files):
    del config, files
    if not NOTEBOOK_OUTPUT_RE.search(markdown):
        return markdown

    page_path = Path(page.file.abs_src_path)
    locale_root = next(
        (parent for parent in page_path.parents if parent.name in {"zh", "en"}),
        None,
    )
    if locale_root is None:
        raise RuntimeError(f"Cannot resolve documentation locale for {page_path}")

    def replace_output(match: re.Match[str]) -> str:
        notebook_path = (locale_root / match.group("notebook")).resolve()
        try:
            notebook_path.relative_to(locale_root.resolve())
        except ValueError as exc:
            raise RuntimeError(
                f"Notebook output path leaves the locale documentation root: {notebook_path}"
            ) from exc
        if not notebook_path.is_file():
            raise RuntimeError(f"Notebook output source does not exist: {notebook_path}")

        output = _notebook_cell_html(notebook_path, match.group("cell_id"))
        wrapped = f'<div class="molop-notebook-output" markdown="0">\n{output}\n</div>'
        indent = match.group("indent")
        return "\n".join(f"{indent}{line}" if line else indent for line in wrapped.splitlines())

    return NOTEBOOK_OUTPUT_RE.sub(replace_output, markdown)


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

        batch_page = Path(config.site_dir) / locale_prefix / "guides" / "batch" / "index.html"
        batch_output = batch_page.read_text(encoding="utf-8")
        if (
            'class="dataframe"' not in batch_output
            or "DiskStorage.FilePath" not in batch_output
            or "Energy.total_energy.hartree" not in batch_output
        ):
            raise RuntimeError(
                f"Complete notebook summary table is missing from page: {batch_page}"
            )
