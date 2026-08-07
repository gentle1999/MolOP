# Documentation and MolOP versions

This site is built from the current MolOP source tree and keeps snapshots for stable release tags. Use the version selector to switch among the latest stable release, previous releases, and the `main` development build. The banner at the top of every page identifies the source version, Git ref, commit, and latest release represented by the current snapshot. Include those values when reporting a documentation or API mismatch.

## Documentation states

| State | Meaning |
| --- | --- |
| Development documentation | The source is an unreleased commit on `main` or a locally modified checkout. It may describe APIs that the latest PyPI release does not provide. |
| Release documentation | The source exactly matches a stable release tag and contains no additional working-tree changes. |

Stable snapshots use `/<version>/` paths, `latest` points to the newest stable release, and `main`/`dev` point to the current development build. The site root defaults to `latest`, so old unversioned root links resolve to the newest stable documentation while versioned links remain stable.

Snapshots are built from stable release tags, which update `latest`. Pre-release tags are outside this versioning system. Tags created before this versioning workflow was enabled do not have snapshots automatically; use the documentation deployment workflow's `workflow_dispatch` to backfill a selected stable tag.

For a backfill, set `source_ref` to the target stable tag, such as `v0.2.2`. Set `update_latest` to `true` only when that tag is the newest stable release; leave it `false` for older tags so they cannot move `latest` backwards.

## Check an installed version

```bash
python -c "import molop; print(molop.__version__)"
```

If that output differs from the source version in the banner, check the latest release shown there:

- For a PyPI installation, use the API contract from the corresponding release tag.
- To verify a new mainline API, install the commit shown in the banner instead of assuming the latest PyPI package already contains it.

```bash
uv add "molop @ git+https://github.com/gentle1999/MolOP.git@COMMIT"
```

Replace `COMMIT` with the commit from the banner. Pin a release version or full commit for production and reproducible research; do not depend on a moving `main` branch.

## Plugin compatibility

Third-party reader/writer plugins should declare the tested MolOP version range in `pyproject.toml`. If a plugin depends on a new interface described only by development documentation, pin the mainline commit while developing. Before publishing the plugin, wait for that interface to enter a MolOP release and update the dependency lower bound.

See [Reader/writer plugin development](developer/plugins.md) for the registration contract.
