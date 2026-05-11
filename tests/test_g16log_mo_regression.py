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
