from pathlib import Path

from molop import AutoParser


FIXTURE = (
    Path(__file__).resolve().parent
    / "test_files"
    / "orca"
    / "output_files"
    / "local"
    / "H2_sp_orca.out"
)


def test_representative_orca_output_is_in_regular_regression() -> None:
    batch = AutoParser(FIXTURE, parser_detection="auto", n_jobs=1, only_last_frame=True)

    assert len(batch) == 1
    assert batch[0].detected_format_id == "orcaout"
    assert batch[0].qm_software == "ORCA"
    assert batch[0].status.normal_terminated is True
