from pathlib import Path

from molop.io.logic.orca.log.locators import locate_orca_job_frames, locate_orca_jobs


ORCA_FIXTURE_ROOT = Path(__file__).resolve().parent / "test_files" / "orca"


def test_orca_optimization_frames_are_exact_source_partition() -> None:
    source = (ORCA_FIXTURE_ROOT / "opt_orca.out").read_text(encoding="utf-8")
    jobs = locate_orca_jobs(source)

    assert len(jobs) == 1
    frames = locate_orca_job_frames(source, jobs[0])
    assert len(frames) == 5
    assert "CARTESIAN COORDINATES (ANGSTROEM)" in frames[0].text(source)
    assert "OPTIMIZATION RUN DONE" in frames[-1].text(source)
    assert "".join(frame.text(source) for frame in frames) == source


def test_orca_job_locator_does_not_search_for_repeated_frame_text() -> None:
    coordinate_block = (
        "------------------\nCARTESIAN COORDINATES (ANGSTROEM)\n------------------\nH 0.0 0.0 0.0\n"
    )
    source = (
        "*** JOB NUMBER 1 ***\n"
        f"{coordinate_block}same-result\n"
        f"{coordinate_block}same-result\n"
        "*** JOB NUMBER 2 ***\n"
        f"{coordinate_block}same-result\n"
    )

    jobs = locate_orca_jobs(source)
    assert len(jobs) == 2
    assert "JOB NUMBER 1" in jobs[0].text(source)
    assert jobs[1].text(source).startswith("*** JOB NUMBER 2 ***")

    first_job_frames = locate_orca_job_frames(source, jobs[0])
    assert len(first_job_frames) == 2
    assert "".join(frame.text(source) for frame in first_job_frames) == jobs[0].text(source)
    assert "".join(job.text(source) for job in jobs) == source


def test_orca_single_point_locator_keeps_preamble_with_observation() -> None:
    source = (ORCA_FIXTURE_ROOT / "h2_grad_orca.out").read_text(encoding="utf-8")
    jobs = locate_orca_jobs(source)
    frames = locate_orca_job_frames(source, jobs[0])

    assert len(jobs) == len(frames) == 1
    assert frames[0].start_char == 0
    assert frames[0].end_char == len(source)
    assert "Program Version" in frames[0].text(source)
    assert "CARTESIAN GRADIENT" in frames[0].text(source)


def test_orca_job_locator_keeps_jobs_without_parseable_coordinate_frames() -> None:
    source = (
        "* O   R   C   A *\n"
        "Program Version 6.0.0\n"
        "********** JOB NUMBER 2 **********\n"
        "------------------\n"
        "CARTESIAN COORDINATES (ANGSTROEM)\n"
        "------------------\n"
        "coordinate parsing failed\n"
    )

    jobs = locate_orca_jobs(source)

    assert len(jobs) == 2
    assert locate_orca_job_frames(source, jobs[0]) == ()
    assert locate_orca_job_frames(source, jobs[1]) == ()
    assert "JOB NUMBER 2" in jobs[1].text(source)
