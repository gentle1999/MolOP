from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern


XTB_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][-+]?\d+)?"


class XTBOutputPatterns:
    BANNER = MolOPPattern(content_pattern=r"(?i)\bx\s+T\s+B\b")
    MODERN_VERSION = MolOPPattern(
        content_pattern=(
            r"(?im)^\s*\*\s*xtb\s+version\s+"
            r"(?P<version>\d+(?:\.\d+){1,2}(?:[-+._A-Za-z0-9]*))"
        ),
        content_repeat=0,
    )
    LEGACY_VERSION = MolOPPattern(
        content_pattern=(
            r"(?im)^\s*\|\s*Version\s+"
            r"(?P<version>\d+(?:\.\d+){1,2}(?:\s+[^|]+?)?)\s*\|\s*$"
        ),
        content_repeat=0,
    )
    PROGRAM_CALL = MolOPPattern(
        content_pattern=r"(?im)^\s*program call\s*:\s*(?P<value>.+?)\s*$",
        content_repeat=0,
    )
    COORDINATE_FILE = MolOPPattern(
        content_pattern=r"(?im)^\s*coordinate file\s*:\s*(?P<value>.+?)\s*$",
        content_repeat=0,
    )
    OMP_THREADS = MolOPPattern(
        content_pattern=r"(?im)^\s*omp threads\s*:\s*(?P<value>\d+)\s*$",
        content_repeat=0,
    )
    HAMILTONIAN = MolOPPattern(
        content_pattern=r"(?im)Hamiltonian\s+(?P<value>GFN(?:\d+|[-]?FF)?-xTB)\b",
        content_repeat=0,
    )
    SPACED_HAMILTONIAN = MolOPPattern(
        content_pattern=(
            r"(?im)^\s*\|\s*G\s*F\s*N\s*"
            r"(?:(?P<level>\d+|[-]?F\s*F)\s*)?-\s*x\s*T\s*B\s*\|"
        ),
        content_repeat=0,
    )
    SOLVATION_MODEL = MolOPPattern(
        content_pattern=r"(?im)^\s*\*?\s*Solvation model\s*:\s*(?P<value>\S+)",
        content_repeat=0,
    )
    SOLVENT = MolOPPattern(
        content_pattern=r"(?im)^\s*Solvent\s+(?P<value>\S+)\s*$",
        content_repeat=0,
    )
    ELECTRON_TEMPERATURE = MolOPPattern(
        content_pattern=(
            rf"(?im)^\s*(?:electronic temp\.|T\(el\))\s*:\s*"
            rf"(?P<value>{XTB_FLOAT_PATTERN})"
        ),
        content_repeat=0,
    )
    VIB_TEMPERATURE = MolOPPattern(
        content_pattern=rf"(?im)^\s*(?P<value>{XTB_FLOAT_PATTERN})\s+VIB\b",
        content_repeat=0,
    )
    SOLVENT_TEMPERATURE = MolOPPattern(
        content_pattern=rf"(?im)^\s*Temperature\s+(?P<value>{XTB_FLOAT_PATTERN})\s+K\b",
        content_repeat=0,
    )
    SETUP_CHARGE = MolOPPattern(
        content_pattern=rf"(?im)^\s*charge\s*:\s*(?P<value>{XTB_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    SETUP_SPIN = MolOPPattern(
        content_pattern=rf"(?im)^\s*spin\s*:\s*(?P<value>{XTB_FLOAT_PATTERN})\s*$",
        content_repeat=0,
    )
    TOTAL_CHARGE = MolOPPattern(
        content_pattern=rf"(?im)total charge\s+(?P<value>{XTB_FLOAT_PATTERN})\s+e\b",
        content_repeat=0,
    )
    FINAL_ENERGIES = (
        MolOPPattern(
            content_pattern=(
                rf"(?im)^\s*\|\s*TOTAL ENERGY\s+"
                rf"(?P<value>{XTB_FLOAT_PATTERN})\s+Eh\s*\|"
            ),
            content_repeat=0,
        ),
        MolOPPattern(
            content_pattern=(
                rf"(?im)^\s*::\s*total energy\s+"
                rf"(?P<value>{XTB_FLOAT_PATTERN})\s+Eh\s*::"
            ),
            content_repeat=0,
        ),
        MolOPPattern(
            content_pattern=rf"(?im)^\s*total E\s*:\s*(?P<value>{XTB_FLOAT_PATTERN})\b",
            content_repeat=0,
        ),
        MolOPPattern(
            content_pattern=(
                rf"(?im)^\s*\*\s*total energy\s*:\s*"
                rf"(?P<value>{XTB_FLOAT_PATTERN})\s+Eh"
            ),
            content_repeat=0,
        ),
    )
    TIME_VALUE = (
        r"(?P<days>\d+)\s*d,\s*(?P<hours>\d+)\s*h,\s*"
        rf"(?P<minutes>\d+)\s*min,\s*(?P<seconds>{XTB_FLOAT_PATTERN})\s*sec"
    )
    TOTAL_WALL_TIME = MolOPPattern(
        content_pattern=rf"(?im)^\s*total:\s*$\s*^\s*\*\s*wall-time:\s*{TIME_VALUE}",
        content_repeat=0,
    )
    ANY_WALL_TIME = MolOPPattern(
        content_pattern=rf"(?im)^\s*\*\s*wall-time:\s*{TIME_VALUE}",
        content_repeat=0,
    )
    FINAL_STRUCTURE = MolOPPattern(
        content_pattern=r"(?im)^\s*final structure:\s*$",
        content_repeat=0,
    )
    LEGACY_COORDINATE_HEADER = MolOPPattern(
        content_pattern=r"(?im)^\s*Z\s+AO/shell.*Cart\. coordinates\s*$"
    )
    ENERGY_CONVERGENCE = MolOPPattern(
        content_pattern=(rf"(?im)^\s*energy convergence\s+(?P<value>{XTB_FLOAT_PATTERN})"),
        content_repeat=0,
    )
    ENERGY_CHANGE = MolOPPattern(
        content_pattern=rf"(?im)change\s+(?P<value>{XTB_FLOAT_PATTERN})\s+Eh",
        content_repeat=0,
    )
    GRADIENT_CONVERGENCE = MolOPPattern(
        content_pattern=(rf"(?im)^\s*grad\. convergence\s+(?P<value>{XTB_FLOAT_PATTERN})"),
        content_repeat=0,
    )
    GRADIENT_NORM = MolOPPattern(
        content_pattern=(rf"(?im)^\s*gradient norm\s*:\s*(?P<value>{XTB_FLOAT_PATTERN})\s+Eh"),
        content_repeat=0,
    )
    INDEXED_VALUE = MolOPPattern(
        content_pattern=rf"(?P<index>\d+)\s*:\s*(?P<value>{XTB_FLOAT_PATTERN})",
        content_repeat=0,
    )
    FREQUENCY_HEADER = MolOPPattern(
        content_pattern=r"(?im)^\s*projected vibrational frequencies \(cm.¹\)\s*$",
        content_repeat=0,
    )
    REDUCED_MASS_HEADER = MolOPPattern(content_pattern=r"(?im)^\s*reduced masses \(amu\)\s*$")
    IR_HEADER = MolOPPattern(content_pattern=r"(?im)^\s*IR intensities \(km.mol.¹\)\s*$")
    RAMAN_HEADER = MolOPPattern(content_pattern=r"(?im)^\s*Raman intensities")
    FREQUENCY_ROW = MolOPPattern(
        content_pattern=r"(?im)^\s*eigval\s*:\s*(?P<values>.+)$",
        content_repeat=0,
    )
    FLOAT_TOKEN = MolOPPattern(
        content_pattern=rf"(?P<value>{XTB_FLOAT_PATTERN})",
        content_repeat=0,
    )
    ZPVE = MolOPPattern(
        content_pattern=(
            rf"(?im)^\s*::\s*zero point energy\s+"
            rf"(?P<value>{XTB_FLOAT_PATTERN})\s+Eh"
        ),
        content_repeat=0,
    )
    TOTAL_ENTHALPY = MolOPPattern(
        content_pattern=(
            rf"(?im)^\s*\|\s*TOTAL ENTHALPY\s+"
            rf"(?P<value>{XTB_FLOAT_PATTERN})\s+Eh"
        ),
        content_repeat=0,
    )
    TOTAL_FREE_ENERGY = MolOPPattern(
        content_pattern=(
            rf"(?im)^\s*\|\s*TOTAL FREE ENERGY\s+"
            rf"(?P<value>{XTB_FLOAT_PATTERN})\s+Eh"
        ),
        content_repeat=0,
    )
    THERMO_TOTAL_ROW = MolOPPattern(
        content_pattern=r"(?im)^\s*TOT\s+(?P<values>.+)$",
        content_repeat=0,
    )
    MOLECULAR_MASS = MolOPPattern(
        content_pattern=(rf"(?im)^\s*molecular mass/u\s*:\s*(?P<value>{XTB_FLOAT_PATTERN})"),
        content_repeat=0,
    )
    ROTATION_CONSTANTS = MolOPPattern(
        content_pattern=(
            rf"(?im)^\s*rotational constants/cm.¹\s*:\s*"
            rf"(?P<a>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<b>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<c>{XTB_FLOAT_PATTERN})"
        ),
        content_repeat=0,
    )
    GFN1_CHARGE_HEADER = MolOPPattern(
        content_pattern=(r"(?im)^\s*Mulliken/CM5 charges\s+n\(s\)\s+n\(p\)\s+n\(d\)\s*$"),
        content_repeat=0,
    )
    GFN1_CHARGE_ROW = MolOPPattern(
        content_pattern=(
            rf"^\s*\d+[A-Za-z]{{1,3}}\s+"
            rf"(?P<mulliken>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<cm5>{XTB_FLOAT_PATTERN})"
        )
    )
    GFN2_CHARGE_HEADER = MolOPPattern(
        content_pattern=r"(?im)^\s*#\s+Z\s+(?:sym\s+)?covCN\s+q\s+C6AA",
        content_repeat=0,
    )
    GFN2_CHARGE_ROW = MolOPPattern(
        content_pattern=(
            rf"^\s*\d+\s+\d+\s+[A-Za-z]{{1,3}}\s+{XTB_FLOAT_PATTERN}\s+"
            rf"(?P<charge>{XTB_FLOAT_PATTERN})"
        )
    )
    MODERN_ORBITAL_HEADER = MolOPPattern(
        content_pattern=r"(?im)^\s*\*\s*Orbital Energies and Occupations\s*$",
        content_repeat=0,
    )
    HL_GAP = MolOPPattern(content_pattern=r"(?im)^\s*HL-Gap\b")
    LEGACY_ORBITAL_HEADER = MolOPPattern(
        content_pattern=r"(?im)^\s*eigenvalues\s*$",
        content_repeat=0,
    )
    SCC_ENERGY_HEADER = MolOPPattern(content_pattern=r"(?im)^\s*SCC energy\s*:")
    LEGACY_ORBITAL_ROWS = MolOPPattern(
        content_pattern=(
            r"(?im)^\s*occ\.\s*:\s*(?P<occupancies>.+)$\s*"
            r"^\s*eps\s*:\s*(?P<energies>.+)$"
        ),
        content_repeat=0,
    )
    DIPOLE = MolOPPattern(
        content_pattern=(
            rf"(?im)^\s*full:\s*(?P<x>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<y>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<z>{XTB_FLOAT_PATTERN})"
            rf"(?:\s+(?P<total>{XTB_FLOAT_PATTERN}))?\s*$"
        ),
        content_repeat=0,
    )
    VIP = MolOPPattern(
        content_pattern=(rf"(?im)delta SCC IP \(eV\):\s*(?P<value>{XTB_FLOAT_PATTERN})"),
        content_repeat=0,
    )
    VEA = MolOPPattern(
        content_pattern=(rf"(?im)delta SCC EA \(eV\):\s*(?P<value>{XTB_FLOAT_PATTERN})"),
        content_repeat=0,
    )
    GEI = MolOPPattern(
        content_pattern=(
            rf"(?im)Global electrophilicity index \(eV\):\s*"
            rf"(?P<value>{XTB_FLOAT_PATTERN})"
        ),
        content_repeat=0,
    )
    FUKUI_HEADER = MolOPPattern(
        content_pattern=r"(?im)^\s*#\s+f\(\+\)\s+f\(-\)\s+f\(0\)\s*$",
        content_repeat=0,
    )
    FUKUI_ROW = MolOPPattern(
        content_pattern=(
            rf"^\s*\d+[A-Za-z]{{1,3}}\s+"
            rf"(?P<positive>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<negative>{XTB_FLOAT_PATTERN})\s+"
            rf"(?P<zero>{XTB_FLOAT_PATTERN})"
        )
    )


xtb_output_patterns = XTBOutputPatterns()


__all__ = ["XTB_FLOAT_PATTERN", "XTBOutputPatterns", "xtb_output_patterns"]
