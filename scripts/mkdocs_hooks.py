from __future__ import annotations

from mkdocs import plugins


@plugins.event_priority(-200)
def on_files(files, config):
    _ = config
    for file in files:
        if file.src_uri.endswith(".ipynb"):
            file.is_documentation_page = lambda: True
    return files
