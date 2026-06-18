from __future__ import annotations

import warnings

from mkdocs import plugins
from nbformat.warnings import MissingIDFieldWarning


warnings.filterwarnings("ignore", category=MissingIDFieldWarning)


@plugins.event_priority(-200)
def on_files(files, config):
    _ = config
    for file in files:
        if file.src_uri.startswith("examples/") and file.src_uri.endswith(".ipynb"):
            file.is_documentation_page = lambda: True
    return files
