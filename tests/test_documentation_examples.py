from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from molop import AutoParser
from molop.cli.app import app


EXAMPLE_DIR = Path("docs/assets/examples")
ORCA_EXAMPLE = EXAMPLE_DIR / "water_mp2.out"
USER_DOC_GLOBS = (
    "index.md",
    "getting_started/*.md",
    "guides/*.md",
    "tutorials/*.md",
)
PAIRED_USER_PAGES = (
    "index.md",
    "concepts.md",
    "command_line_interface.md",
    "getting_started/installation.md",
    "getting_started/quickstart.md",
    "getting_started/python-api.md",
    "getting_started/cli.md",
    "guides/parsing.md",
    "guides/results.md",
    "guides/batch.md",
    "guides/filtering.md",
    "guides/conversion.md",
    "guides/structure-recovery.md",
    "guides/cli-recipes.md",
    "guides/troubleshooting.md",
    "tutorials/index.md",
    "tutorials/qm-summary.md",
    "tutorials/energy-csv.md",
    "tutorials/select-results.md",
    "tutorials/export-inputs.md",
    "reference/format_support.md",
    "reference/model_fields.md",
)
NOTEBOOK_NAMES = (
    "01-gaussian-parse-and-inspect.ipynb",
    "02-batch-summary-filter-select.ipynb",
    "03-transform-and-export.ipynb",
)


def test_quickstart_orca_example_exposes_structure_energy_and_summary() -> None:
    batch = AutoParser(ORCA_EXAMPLE, n_jobs=1)
    parsed_file = batch[0]
    frame = parsed_file[-1]

    assert parsed_file.detected_format_id == "orcaout"
    assert len(frame.atoms) == 3
    assert frame.coords.shape == (3, 3)
    assert frame.energies is not None
    assert frame.energies.total_energy is not None
    assert frame.energies.total_energy.m_as("hartree") == pytest.approx(-74.999374598107)
    assert frame.energies.total_energy.m_as("eV") == pytest.approx(-2040.8369503966658)

    summary = batch.to_summary_df(brief=False, flatten_columns=True)
    assert len(summary) == 1
    assert summary.loc[0, "Status.IsNormal"]
    assert summary.loc[0, "Energy.total_energy.hartree"] == pytest.approx(-74.999374598107)


def test_documented_python_batch_operations_are_runnable() -> None:
    batch = AutoParser([ORCA_EXAMPLE], n_jobs=1)

    assert len(batch.filter_state("normal")) == 1
    assert len(batch.filter_by_codec_id("orcaout")) == 1
    rendered = batch.format_transform("xyz", write_to_disk=False)
    xyz = rendered[str(ORCA_EXAMPLE.resolve())]
    assert xyz.splitlines() == [
        "3",
        "comment charge 0 multiplicity 1",
        "O               1.7849140000      1.2624220000      0.5119850000",
        "H               2.6482370000      1.0729290000      0.1316310000",
        "H               1.1831680000      1.2568160000     -0.2388350000",
    ]

    frame = batch[0][-1]
    assert frame.rdmol is not None
    assert frame.rdmol.GetNumAtoms() == 3
    assert frame.rdmol.GetNumBonds() == 2
    assert frame.smiles == "[H]O[H]"
    assert frame.topology_reconstruction_status == "succeeded"

    cml = frame.format_transform("cml")
    assert "<cml" in cml

    gjf = frame.format_transform(
        "gjf",
        link0_commands={"nprocshared": "16", "mem": "32GB"},
        route_section="#p wb97xd/def2tzvp opt freq",
        title_card="optimization and frequency",
        coords_type="cartesian",
        chk=True,
    )
    assert gjf.startswith("%nprocshared=16\n%mem=32GB\n%chk=")
    assert "#p wb97xd/def2tzvp opt freq" in gjf

    orca_input = frame.format_transform(
        "orcainp",
        keywords=("wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP"),
        nprocs=16,
        maxcore=4000,
        blocks={"scf": {"MaxIter": 300, "STABPerform": True}},
    )
    assert orca_input.startswith("! wB97M-V def2-TZVPP")
    assert "%pal\n  nprocs 16\nend" in orca_input
    assert "%maxcore 4000" in orca_input


def test_documented_cli_summary_command_is_runnable(tmp_path: Path) -> None:
    output_path = tmp_path / "summary.csv"
    result = CliRunner().invoke(
        app,
        [
            "-q",
            "parse",
            str(ORCA_EXAMPLE),
            "--n-jobs",
            "1",
            "to-summary-df",
            "--full",
            "--out",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_path.is_file()
    assert "Energy.total_energy.hartree" in output_path.read_text(encoding="utf-8")


def test_documentation_example_provenance_and_license_are_present() -> None:
    source_note = (EXAMPLE_DIR / "SOURCE.txt").read_text(encoding="utf-8")
    license_text = (EXAMPLE_DIR / "LICENSE.cclib").read_text(encoding="utf-8")

    assert "cclib" in source_note
    assert "water_mp2.out" in source_note
    assert "BSD 3-Clause License" in license_text


@pytest.mark.parametrize("locale", ("zh", "en"))
@pytest.mark.parametrize("notebook_name", NOTEBOOK_NAMES)
def test_notebook_metadata_supports_mkdocs_rendering(locale: str, notebook_name: str) -> None:
    notebook = json.loads(
        (Path("docs") / locale / "examples" / notebook_name).read_text(encoding="utf-8")
    )

    assert notebook["metadata"]["kernelspec"] == {
        "display_name": "molop",
        "language": "python",
        "name": "python3",
    }
    assert notebook["metadata"]["language_info"]["name"] == "python"


@pytest.mark.parametrize("relative_path", PAIRED_USER_PAGES)
def test_core_user_documentation_is_bilingual(relative_path: str) -> None:
    assert (Path("docs/zh") / relative_path).is_file()
    assert (Path("docs/en") / relative_path).is_file()


def test_user_main_path_does_not_depend_on_source_checkout() -> None:
    forbidden = ("uv run molop", "tests/test_files", "pyproject.toml", "TQDM_DISABLE")

    for locale in ("zh", "en"):
        locale_root = Path("docs") / locale
        pages = [path for pattern in USER_DOC_GLOBS for path in locale_root.glob(pattern)]
        for page in pages:
            text = page.read_text(encoding="utf-8")
            assert not any(token in text for token in forbidden), page


@pytest.mark.parametrize(
    "relative_path",
    (
        "getting_started/quickstart.md",
        "getting_started/python-api.md",
        "guides/parsing.md",
        "guides/results.md",
        "guides/batch.md",
        "guides/filtering.md",
        "guides/conversion.md",
        "guides/structure-recovery.md",
        "guides/cli-recipes.md",
    ),
)
def test_core_examples_include_output_demonstrations(relative_path: str) -> None:
    for locale, output_markers in (("zh", ("输出", "结果")), ("en", ("output", "result"))):
        text = (Path("docs") / locale / relative_path).read_text(encoding="utf-8")
        assert "```" in text
        searchable_text = text if locale == "zh" else text.lower()
        assert any(marker in searchable_text for marker in output_markers)
