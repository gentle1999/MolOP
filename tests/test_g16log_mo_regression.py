from pathlib import Path

from molop import AutoParser


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "Lnew3_MNNM_final.log"


def test_default_g16log_parser_preserves_open_shell_mo_counts() -> None:
    batch = AutoParser(str(FIXTURE), n_jobs=1)
    frame = batch[0][-1]
    mos = frame.molecular_orbitals

    assert mos is not None
    assert len(mos.alpha_energies) == 1484
    assert len(mos.beta_energies) == 1484
    assert len(mos.alpha_occupancies) == 1484
    assert len(mos.beta_occupancies) == 1484
    assert mos.HOMO_id == 340
    assert mos.beta_HOMO_id == 338
    assert len(mos.SOMO_ids) == 2


MERGED_VALUE_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "6-INT4-Opt.log"
MISSING_SEPARATOR_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "119907.log"


def test_default_g16log_parser_handles_concatenated_orbital_energies() -> None:
    batch = AutoParser(str(MERGED_VALUE_FIXTURE), n_jobs=1)
    frame = batch[0][-1]
    mos = frame.molecular_orbitals

    assert mos is not None
    assert len(mos.alpha_energies) == len(mos.alpha_occupancies)
    assert len(mos.alpha_energies) > 1000
    assert mos.HOMO_id is not None
    assert mos.electronic_state == "1-A"
    assert len(mos.beta_energies) == len(mos.beta_occupancies)


def test_g16log_orbital_energies_without_separators_remain_aligned() -> None:
    file_model = AutoParser(
        str(MISSING_SEPARATOR_FIXTURE),
        parser_detection="g16log",
        n_jobs=1,
    )[0]
    mos = file_model[0].molecular_orbitals

    assert mos is not None
    assert mos.alpha_energies[:5].to("hartree").magnitude.tolist() == [
        -255.98335,
        -255.98334,
        -255.98328,
        -29.93880,
        -29.93880,
    ]
    assert len(mos.alpha_energies) == 429
    assert len(mos.alpha_occupancies) == 429
    assert len(mos.alpha_symmetries) == 429
    assert len(mos.beta_energies) == 429
    assert len(mos.beta_occupancies) == 429
    assert len(mos.beta_symmetries) == 429
    assert mos.HOMO_id == 123
    assert mos.beta_HOMO_id == 122
    # Initial-guess orbital 4 is ``(E)``; the converged population-analysis
    # block reports it as ``(?A)``. Only the latter belongs in the frame model.
    assert mos.alpha_symmetries[3] == "(?A)"
    assert mos.beta_symmetries[3] == "(?A)"
    assert "(?B)" in mos.alpha_symmetries
    assert "(?B)" in mos.beta_symmetries
