from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from click.testing import CliRunner
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from molop import AutoParser, molopconfig
from molop.cli.app import app


EXAMPLE_DIR = Path("docs/assets/examples")
ORCA_EXAMPLE = EXAMPLE_DIR / "water_mp2.out"
METAL_COMPLEX_EXAMPLE = EXAMPLE_DIR / "mn_complex_sp.log"
USER_DOC_GLOBS = (
    "index.md",
    "getting_started/*.md",
    "guides/*.md",
    "tutorials/*.md",
    "reference/*.md",
    "reference/api/*.md",
    "reference/behavior/*.md",
    "reference/formats/*.md",
)
PAIRED_USER_PAGES = (
    "index.md",
    "getting_started/concepts.md",
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
    "tutorials/transition-states.md",
    "tutorials/export-inputs.md",
    "reference/cli.md",
    "reference/behavior/serialization.md",
    "reference/format_support.md",
    "reference/model_fields.md",
    "developer/contracts/api.md",
)
NOTEBOOK_NAMES = (
    "01-gaussian-parse-and-inspect.ipynb",
    "02-batch-summary-filter-select.ipynb",
    "03-transform-and-export.ipynb",
    "04-transition-state-endpoints.ipynb",
    "05-substituent-replacement.ipynb",
)
GENERATED_VISUALIZATION_ASSETS = (
    "mn_complex_graph_reconstruction.svg",
    "ts_imaginary_mode.svg",
    "ts_endpoints_difference.svg",
    "ts_endpoints_additional_sampling.svg",
)

OUTPUT_SECTION_RE = re.compile(
    r"^#{1,6}\s+(?:expected output|sample output|output|outputs|输出|预期输出|示例输出|运行结果)\s*:?[：:]?\s*$",
    re.IGNORECASE,
)
NAKED_OUTPUT_FENCE_RE = re.compile(
    r"^(?:expected output|sample output|output|outputs|输出|预期输出|示例输出|运行结果)\s*[:：]?\s*\n\s*```",
    re.IGNORECASE | re.MULTILINE,
)
PYPI_RELATIVE_LINK_RE = re.compile(r"\]\((?!https?://|mailto:|#)[^)]+\)")
OUTPUT_FENCE_LANGUAGES = {"text", "console", "csv", "tree"}
EXECUTABLE_FENCE_RE = re.compile(r"^\s*```(python|bash|console)\s*$")
OUTPUT_PRODUCING_CODE_RE = re.compile(
    r"(?:\bprint\s*\(|\bdisplay\s*\(|\.to_csv\s*\(|"
    r"\bformat_transform\s*\(|^\s*molop(?:\s|$))",
    re.MULTILINE,
)


def _markdown_lines_without_fenced_content(text: str) -> list[str | None]:
    """Keep headings visible while ignoring illustrative Markdown snippets."""
    result: list[str | None] = []
    in_fence = False
    for line in text.splitlines():
        if line.lstrip().startswith("```") or line.lstrip().startswith("````"):
            in_fence = not in_fence
            result.append(None)
        else:
            result.append(None if in_fence else line)
    return result


def _executable_markdown_blocks(text: str) -> list[tuple[int, int, str]]:
    lines = text.splitlines()
    blocks: list[tuple[int, int, str]] = []
    index = 0
    while index < len(lines):
        match = EXECUTABLE_FENCE_RE.match(lines[index])
        if not match:
            index += 1
            continue

        start = index
        index += 1
        content: list[str] = []
        while index < len(lines) and not lines[index].lstrip().startswith("```"):
            content.append(lines[index])
            index += 1
        blocks.append((start, index, "\n".join(content)))
        index += 1
    return blocks


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
    assert summary.shape == (1, 23)
    assert summary.loc[0, "Status.IsNormal"]
    assert summary.loc[0, "Energy.total_energy.hartree"] == pytest.approx(-74.999374598107)


def test_mixed_qm_outputs_share_one_batch_container() -> None:
    batch = AutoParser([METAL_COMPLEX_EXAMPLE, ORCA_EXAMPLE], n_jobs=1)
    rows = {parsed_file.detected_format_id: parsed_file[-1] for parsed_file in batch}

    assert set(rows) == {"g16log", "orcaout"}
    assert (rows["g16log"].qm_software, rows["g16log"].method) == ("Gaussian", "DFT")
    assert (rows["orcaout"].qm_software, rows["orcaout"].method) == ("ORCA", "MP2")
    assert rows["g16log"].energies.total_energy.m_as("hartree") == pytest.approx(-2182.472195)
    assert rows["orcaout"].energies.total_energy.m_as("hartree") == pytest.approx(-74.999374598107)
    assert batch.to_summary_df(brief=False, flatten_columns=True).shape == (2, 23)


def test_homepage_leads_with_program_independent_processing() -> None:
    expected_content = {
        "zh": (
            "MolOP 是面向计算化学文件的 Python 库和命令行工具",
            "统一整理为 `batch -> file -> frame` 对象模型",
            "## 数据流概览",
        ),
        "en": (
            "MolOP is a Python library and command-line tool for computational chemistry files",
            "maps them into a common `batch -> file -> frame` object model",
            "## Data-flow overview",
        ),
    }
    for locale, expected in expected_content.items():
        text = (Path("docs") / locale / "index.md").read_text(encoding="utf-8")
        normalized_text = " ".join(text.split())
        first_h2 = next(line for line in text.splitlines() if line.startswith("## "))
        assert first_h2 == expected[2]
        assert all(" ".join(fragment.split()) in normalized_text for fragment in expected)
        assert "Gaussian 16" not in text
        assert "mn_complex_graph_reconstruction.svg" not in text

    readmes = {
        Path("README.zh.md"): ("## 统一读取不同量化软件文件", "## 从坐标恢复分子图"),
        Path("README.md"): (
            "## Read different quantum-chemistry files through one model",
            "## Recover molecular graphs from coordinates",
        ),
    }
    for path, (core_heading, graph_heading) in readmes.items():
        text = path.read_text(encoding="utf-8")
        assert core_heading in text
        assert graph_heading not in text


def test_homepage_mermaid_maps_inputs_to_one_batch_and_downstream_tasks() -> None:
    localized_labels = {
        "zh": ("不同来源与文件格式", "统一的下游处理"),
        "en": ("Different sources and file formats", "Common downstream processing"),
    }

    for locale, labels in localized_labels.items():
        text = (Path("docs") / locale / "index.md").read_text(encoding="utf-8")
        assert text.count("```mermaid") == 1
        assert "flowchart TB" in text
        assert "MolOP AutoParser" in text
        assert "FileBatchModelDisk" in text
        assert "future --> parser" in text
        assert "parser --> batch" in text
        assert all(
            f"batch --> {task}" in text for task in ("results", "summary", "export", "workflow")
        )
        assert all(label in text for label in labels)

    mkdocs_config = Path("mkdocs.yml").read_text(encoding="utf-8")
    assert "- name: mermaid" in mkdocs_config
    assert "https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.min.js" in mkdocs_config


def test_user_installation_uses_the_pypi_package() -> None:
    installation_pages = (
        Path("README.md"),
        Path("README.zh.md"),
        Path("docs/en/getting_started/installation.md"),
        Path("docs/zh/getting_started/installation.md"),
    )
    all_user_docs = [Path("README.md"), Path("README.zh.md"), *Path("docs").rglob("*.md")]

    for page in installation_pages:
        assert "pip install molop" in page.read_text(encoding="utf-8")

    forbidden_release_text = (
        'python -m pip install "git+https://github.com/gentle1999/MolOP.git"',
        "MolOP is not currently published on PyPI",
        "MolOP 当前尚未发布到 PyPI",
    )
    for page in all_user_docs:
        text = page.read_text(encoding="utf-8")
        assert not any(fragment in text for fragment in forbidden_release_text), page


def test_readmes_expose_synchronized_verifiable_status_badges() -> None:
    badge_sources = (
        "img.shields.io/pypi/v/molop.svg",
        "img.shields.io/pypi/pyversions/molop.svg",
        "img.shields.io/badge/typing-typed-blue.svg",
        "img.shields.io/pypi/status/molop.svg",
        "img.shields.io/pypi/wheel/molop.svg",
        "img.shields.io/pypi/dm/molop.svg",
        "actions/workflows/ci.yaml/badge.svg?branch=main",
        "actions/workflows/docs-deploy.yml/badge.svg?branch=main",
        "img.shields.io/github/license/gentle1999/MolOP.svg",
        "img.shields.io/github/last-commit/gentle1999/MolOP.svg",
        "img.shields.io/github/issues/gentle1999/MolOP.svg",
        "img.shields.io/github/stars/gentle1999/MolOP.svg",
        "img.shields.io/github/forks/gentle1999/MolOP.svg",
    )

    badge_blocks: list[list[str]] = []
    for readme in (Path("README.md"), Path("README.zh.md")):
        text = readme.read_text(encoding="utf-8")
        assert all(source in text for source in badge_sources)
        assert "codecov" not in text.lower()
        lines = text.splitlines()
        badge_blocks.append([line for line in lines if line.startswith("[![")])

    assert badge_blocks[0] == badge_blocks[1]


def test_readme_links_are_portable_to_pypi() -> None:
    for readme in (Path("README.md"), Path("README.zh.md")):
        text = readme.read_text(encoding="utf-8")
        assert PYPI_RELATIVE_LINK_RE.search(text) is None, readme


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


def test_documented_batch_grid_image_is_runnable() -> None:
    batch = AutoParser([ORCA_EXAMPLE], n_jobs=1)

    grid = batch.draw_grid_image(maxMols=16, n_jobs=1)

    assert isinstance(grid, str)
    assert "<svg" in grid[:200]


def test_batch_grid_image_uses_dof_drawer_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    import rdkit_dof

    calls: list[dict[str, object]] = []

    def fake_dof_drawer(molecules, **kwargs):
        calls.append({"molecules": molecules, **kwargs})
        return "<svg data-drawer='rdkit-dof'></svg>"

    monkeypatch.setattr(rdkit_dof, "MolsToGridDofImage", fake_dof_drawer)
    batch = AutoParser([ORCA_EXAMPLE], n_jobs=1)

    grid = batch.draw_grid_image(maxMols=16, n_jobs=1)

    assert grid == "<svg data-drawer='rdkit-dof'></svg>"
    assert calls and calls[0]["use_svg"] is True
    assert calls[0]["return_image"] is False


def test_metal_complex_reconstruction_exposes_molgr_advantage(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(molopconfig, "graph_reconstruction_backend", "cpp")
    monkeypatch.setattr(molopconfig, "make_dative_bonds", True)

    frame = AutoParser(METAL_COMPLEX_EXAMPLE, n_jobs=1)[0][-1]
    raw_xyz = Chem.MolFromXYZBlock(frame.to_XYZ())
    assert raw_xyz is not None
    assert raw_xyz.GetNumAtoms() == 32
    assert raw_xyz.GetNumBonds() == 0

    distance_connectivity = Chem.Mol(raw_xyz)
    rdDetermineBonds.DetermineConnectivity(distance_connectivity)
    distance_mn = next(
        atom for atom in distance_connectivity.GetAtoms() if atom.GetSymbol() == "Mn"
    )
    distance_mn_bonds = list(distance_mn.GetBonds())
    assert distance_connectivity.GetNumBonds() == 36
    assert len(distance_mn_bonds) == 8
    assert all(bond.GetBondType() == Chem.BondType.SINGLE for bond in distance_mn_bonds)

    with pytest.raises(ValueError, match="atomic number 25 has no valences defined"):
        determine_bonds = Chem.Mol(raw_xyz)
        rdDetermineBonds.DetermineBonds(determine_bonds, charge=frame.charge)

    reconstructed = frame.rdmol
    assert reconstructed is not None
    assert frame.formula == "C12H15MnO3P+"
    assert reconstructed.GetNumAtoms() == 32
    assert reconstructed.GetNumBonds() == 36
    assert frame.topology_reconstruction_backend == "cpp"
    assert frame.topology_reconstruction_status == "succeeded"

    reconstructed_mn = next(atom for atom in reconstructed.GetAtoms() if atom.GetSymbol() == "Mn")
    reconstructed_mn_bonds = list(reconstructed_mn.GetBonds())
    assert reconstructed_mn.GetFormalCharge() == 1
    assert reconstructed_mn.GetDegree() == 8
    assert len(reconstructed_mn_bonds) == 8
    assert all(bond.GetBondType() == Chem.BondType.DATIVE for bond in reconstructed_mn_bonds)


def test_documented_visualization_assets_are_generated_svg() -> None:
    source_notes = Path("docs/assets/examples/SOURCE.txt").read_text(encoding="utf-8")
    documentation = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (
            Path("docs/zh/guides/structure-recovery.md"),
            Path("docs/en/guides/structure-recovery.md"),
            Path("docs/zh/tutorials/transition-states.md"),
            Path("docs/en/tutorials/transition-states.md"),
        )
    )

    for filename in GENERATED_VISUALIZATION_ASSETS:
        svg = (Path("docs/assets/examples") / filename).read_text(encoding="utf-8")
        assert "<svg" in svg[:500]
        assert filename in source_notes
        assert filename in documentation

    assert "scripts/generate_doc_images.py" in source_notes


@pytest.mark.parametrize("locale", ("zh", "en"))
@pytest.mark.parametrize(
    "notebook_name",
    ("01-gaussian-parse-and-inspect.ipynb", "02-batch-summary-filter-select.ipynb"),
)
def test_notebooks_show_the_complete_summary_table(locale: str, notebook_name: str) -> None:
    notebook = json.loads(
        (Path("docs") / locale / "examples" / notebook_name).read_text(encoding="utf-8")
    )
    sources = ["".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"]
    html_outputs = [
        "".join(output["data"]["text/html"])
        for cell in notebook["cells"]
        for output in cell.get("outputs", [])
        if "text/html" in output.get("data", {})
    ]

    summary_sources = [source for source in sources if "summary = batch.to_summary_df" in source]
    assert summary_sources
    assert any("display(summary)" in source for source in summary_sources)
    assert all("summary[[" not in source for source in summary_sources)
    assert any(
        "<table" in html
        and "<th>DiskStorage.FilePath</th>" in html
        and "<th>Energy.total_energy.hartree</th>" in html
        for html in html_outputs
    )


def test_summary_entry_points_do_not_replace_the_complete_table_with_a_column_slice() -> None:
    pages = (
        Path("README.md"),
        Path("README.zh.md"),
        Path("docs/en/getting_started/quickstart.md"),
        Path("docs/zh/getting_started/quickstart.md"),
        Path("docs/en/reference/api/filebatchmodeldisk.md"),
        Path("docs/zh/reference/api/filebatchmodeldisk.md"),
    )

    for page in pages:
        text = page.read_text(encoding="utf-8")
        assert "summary[[" not in text
        assert "02-batch-summary-filter-select" in text


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
    assert "mn_complex_sp.log" in source_note
    assert "MnCO3C6H6PMe3-mod2-sp-smd-DSDPBEP86d3.log" in source_note
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


@pytest.mark.parametrize("locale", ("zh", "en"))
@pytest.mark.parametrize("notebook_name", NOTEBOOK_NAMES)
def test_notebook_code_cells_are_preexecuted_with_outputs(locale: str, notebook_name: str) -> None:
    notebook = json.loads(
        (Path("docs") / locale / "examples" / notebook_name).read_text(encoding="utf-8")
    )
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]

    assert code_cells
    for cell in code_cells:
        assert isinstance(cell["execution_count"], int)
        assert cell["outputs"]
        assert not any(output["output_type"] == "error" for output in cell["outputs"])


@pytest.mark.parametrize("relative_path", PAIRED_USER_PAGES)
def test_core_user_documentation_is_bilingual(relative_path: str) -> None:
    assert (Path("docs/zh") / relative_path).is_file()
    assert (Path("docs/en") / relative_path).is_file()


def test_all_markdown_pages_have_bilingual_counterparts() -> None:
    english = {path.relative_to(Path("docs/en")) for path in Path("docs/en").rglob("*.md")}
    chinese = {path.relative_to(Path("docs/zh")) for path in Path("docs/zh").rglob("*.md")}

    assert english == chinese


def test_bilingual_markdown_pages_keep_the_same_structural_depth() -> None:
    for english_path in Path("docs/en").rglob("*.md"):
        relative_path = english_path.relative_to(Path("docs/en"))
        chinese_path = Path("docs/zh") / relative_path
        english_text = english_path.read_text(encoding="utf-8")
        chinese_text = chinese_path.read_text(encoding="utf-8")

        english_heading_levels = [
            len(match.group(1))
            for match in re.finditer(r"^(#{1,6})\s+", english_text, re.MULTILINE)
        ]
        chinese_heading_levels = [
            len(match.group(1))
            for match in re.finditer(r"^(#{1,6})\s+", chinese_text, re.MULTILINE)
        ]

        assert english_heading_levels == chinese_heading_levels, relative_path
        assert english_text.count("```") == chinese_text.count("```"), relative_path
        assert english_text.count("??? ") == chinese_text.count("??? "), relative_path


def test_user_workflow_pages_include_folded_or_notebook_rendered_results() -> None:
    for locale in ("zh", "en"):
        locale_root = Path("docs") / locale
        pages = [path for pattern in USER_DOC_GLOBS for path in locale_root.glob(pattern)]
        for page in pages:
            text = page.read_text(encoding="utf-8")
            lines = text.splitlines()
            for start, end, content in _executable_markdown_blocks(text):
                if not OUTPUT_PRODUCING_CODE_RE.search(content):
                    continue
                assert any(
                    line.strip().startswith(("??? example", "<!-- notebook-output:"))
                    for line in lines[end + 1 : end + 21]
                ), f"Executable output is neither folded nor notebook-rendered: {page}:{start + 1}"


def test_writer_reference_examples_are_self_contained() -> None:
    writer_pages = ("cml.md", "fakeg.md", "gjf.md", "orcainp.md")

    for locale in ("zh", "en"):
        for page_name in writer_pages:
            page = Path("docs") / locale / "reference/formats" / page_name
            text = page.read_text(encoding="utf-8")
            assert "from molop import AutoParser" in text, page
            assert "format_transform(" in text, page
            assert "??? example" in text, page


def test_user_documentation_avoids_patronizing_language() -> None:
    disallowed_terms = (
        "小白",
        "傻瓜式",
        "无脑",
        "闭眼照做",
        "零基础也能",
        "idiot-proof",
        "beginner-proof",
        "even a beginner",
        "just copy and paste",
        "simply copy and paste",
    )

    for page in [Path("README.md"), Path("README.zh.md"), *Path("docs").rglob("*.md")]:
        text = page.read_text(encoding="utf-8").lower()
        assert not any(term.lower() in text for term in disallowed_terms), page


def test_api_contracts_stay_task_oriented() -> None:
    expected_sections = {
        "en": (
            "## At a glance",
            "## `AutoParser`",
            "## Frame selector",
            "## `format_transform`",
            "## `to_summary_df`",
            "## `parallel_execute`",
            "## CLI and serialization",
            "## Errors and boundaries",
        ),
        "zh": (
            "## 快速索引",
            "## `AutoParser`",
            "## Frame 选择",
            "## `format_transform`",
            "## `to_summary_df`",
            "## `parallel_execute`",
            "## CLI 与序列化",
            "## 错误与边界",
        ),
    }
    parser_implementation_terms = (
        "_quick_check_file_format",
        "_locate_segments",
        "LocatedSourceSegment",
        "SourceSpan",
    )

    for locale, sections in expected_sections.items():
        text = (Path("docs") / locale / "developer/contracts/api.md").read_text(encoding="utf-8")
        positions = [text.index(section) for section in sections]
        assert positions == sorted(positions)
        assert not any(term in text for term in parser_implementation_terms)


def test_markdown_output_sections_are_collapsed() -> None:
    for path in [Path("README.md"), Path("README.zh.md"), *Path("docs").rglob("*.md")]:
        text = path.read_text(encoding="utf-8")
        visible_lines = _markdown_lines_without_fenced_content(text)
        for index, line in enumerate(visible_lines):
            if line is None or not OUTPUT_SECTION_RE.fullmatch(line.strip()):
                continue

            section_level = len(line) - len(line.lstrip("#"))
            section_lines: list[str] = []
            for candidate in visible_lines[index + 1 :]:
                if candidate is not None and re.match(r"^#{1,6}\s+", candidate):
                    candidate_level = len(candidate) - len(candidate.lstrip("#"))
                    if candidate_level <= section_level:
                        break
                if candidate is not None:
                    section_lines.append(candidate)

            assert any(
                re.match(r"^\?\?\?\s+example\b", candidate) for candidate in section_lines
            ), f"Output section is not collapsed: {path}:{index + 1}"


def test_markdown_has_no_naked_output_fences() -> None:
    for path in [Path("README.md"), Path("README.zh.md"), *Path("docs").rglob("*.md")]:
        text = path.read_text(encoding="utf-8")
        assert NAKED_OUTPUT_FENCE_RE.search(text) is None, path


def test_readme_output_fences_are_inside_details() -> None:
    for path in (Path("README.md"), Path("README.zh.md")):
        details_depth = 0
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if "<details>" in line:
                details_depth += 1
            fence_match = re.match(r"^\s*```([A-Za-z0-9_-]*)\s*$", line)
            if fence_match and fence_match.group(1).lower() in OUTPUT_FENCE_LANGUAGES:
                assert details_depth > 0, f"README output is not collapsed: {path}:{line_number}"
            if "</details>" in line:
                details_depth -= 1

        assert details_depth == 0


def test_notebooks_render_saved_outputs_after_ci_execution() -> None:
    mkdocs_config = Path("mkdocs.yml").read_text(encoding="utf-8")
    assert "- mkdocs-jupyter:" in mkdocs_config
    assert "execute: false" in mkdocs_config

    for workflow in (Path(".github/workflows/docs-deploy.yml"),):
        workflow_text = workflow.read_text(encoding="utf-8")
        assert "nbconvert" in workflow_text
        assert "--execute" in workflow_text
        assert "--inplace" in workflow_text


def test_documentation_exposes_source_and_release_versions() -> None:
    mkdocs_config = Path("mkdocs.yml").read_text(encoding="utf-8")
    assert "custom_dir: overrides" in mkdocs_config
    assert "- 文档版本: developer/release/versioning.md" in mkdocs_config
    assert "provider: mike" in mkdocs_config

    template = Path("overrides/main.html").read_text(encoding="utf-8")
    assert 'class="molop-docs-version"' in template
    assert "source_version" in template
    assert "release_version" in template

    hook = Path("scripts/mkdocs_hooks.py").read_text(encoding="utf-8")
    for variable in (
        "MOLOP_DOCS_SOURCE_VERSION",
        "MOLOP_DOCS_RELEASE_VERSION",
        "MOLOP_DOCS_COMMIT",
        "MOLOP_DOCS_REF",
        "MOLOP_DOCS_CHANNEL",
    ):
        assert variable in hook
    assert "MIKE_DOCS_VERSION" in hook

    for locale in ("zh", "en"):
        version_page = Path("docs") / locale / "developer/release/versioning.md"
        assert version_page.is_file()
        assert "commit" in version_page.read_text(encoding="utf-8").lower()

    for workflow in (Path(".github/workflows/docs-deploy.yml"),):
        workflow_text = workflow.read_text(encoding="utf-8")
        assert "fetch-depth: 0" in workflow_text

    deploy_workflow = Path(".github/workflows/docs-deploy.yml").read_text(encoding="utf-8")
    assert 'tags:\n      - "v*"' in deploy_workflow
    assert "mike deploy" in deploy_workflow
    assert "main dev" in deploy_workflow
    assert '"${DOCS_VERSION}" latest' in deploy_workflow
    assert "mike set-default --push latest" in deploy_workflow
    assert "ref: ${{ inputs.source_ref || github.sha }}" in deploy_workflow
    assert 'echo "MOLOP_DOCS_COMMIT=$commit" >> "$GITHUB_ENV"' in deploy_workflow
    assert "MOLOP_DOCS_COMMIT: ${{ steps.docs-source.outputs.commit }}" in deploy_workflow
    assert "git archive gh-pages" in deploy_workflow
    assert "git archive origin/gh-pages" not in deploy_workflow
    assert "update_latest" in deploy_workflow
    assert "--with mike==2.2.0" in deploy_workflow
    assert "git show origin/main:.github/pages-root-404.html" in deploy_workflow
    assert deploy_workflow.index("name: Resolve documentation version") < deploy_workflow.index(
        "name: Build docs"
    )
    assert "MOLOP_DOCS_SOURCE_VERSION: ${{ steps.docs-version.outputs.source_version }}" in (
        deploy_workflow
    )
    assert "MOLOP_DOCS_CHANNEL: ${{ steps.docs-version.outputs.channel }}" in deploy_workflow

    root_404 = Path(".github/pages-root-404.html").read_text(encoding="utf-8")
    assert "latest/" in root_404
    assert "window.location.replace" in root_404


def test_batch_guides_embed_the_complete_notebook_table() -> None:
    directive = (
        "<!-- notebook-output: examples/02-batch-summary-filter-select.ipynb#batch-summary -->"
    )

    for locale in ("zh", "en"):
        guide = (Path("docs") / locale / "guides" / "batch.md").read_text(encoding="utf-8")
        assert directive in guide

    hook = Path("scripts/mkdocs_hooks.py").read_text(encoding="utf-8")
    assert "_notebook_cell_html" in hook
    assert "Complete notebook summary table is missing" in hook


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
        "tutorials/transition-states.md",
    ),
)
def test_core_examples_include_output_demonstrations(relative_path: str) -> None:
    for locale, output_markers in (("zh", ("输出", "结果")), ("en", ("output", "result"))):
        text = (Path("docs") / locale / relative_path).read_text(encoding="utf-8")
        assert "```" in text
        searchable_text = text if locale == "zh" else text.lower()
        assert any(marker in searchable_text for marker in output_markers)
