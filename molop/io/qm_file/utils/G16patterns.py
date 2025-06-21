"""
Author: TMJ
Date: 2025-02-16 19:20:04
LastEditors: TMJ
LastEditTime: 2025-02-17 12:41:34
Description: 请填写简介
"""


from molop.io.bases.BasePattern import MolOPPattern


class G16InputPatterns:
    """
    G16 input file patterns
    """


class G16LogPatterns:
    """
    G16 log file patterns
    """

    VERSION = MolOPPattern(
        start_pattern=r"^\s\*+\n",
        end_pattern=r"^\s\*+\n",
        content_pattern=r"(Gaussian\s+\d+:\s*[^\s]*\s*\d+-[A-Za-z]{3}-\d{4}\n\s+\d+-[A-Za-z]{3}-\d{4})",
        description="The exact version of Gaussian used. link 1",
    )
    OPTIONS = MolOPPattern(
        start_pattern=r"^\s*\*+\n",
        end_pattern=r"^\s*\-+\n",
        description="The options used in the Gaussian calculation. Such as: MEMORY, CPU, etc. link 1",
    )
    KEYWORDS = MolOPPattern(
        start_pattern=r"^\s*\-+\n",
        end_pattern=r"^\s*\-+\n",
        description="The keywords used in the Gaussian calculation. link 1",
    )
    TITLE = MolOPPattern(
        start_pattern=r"^\s*\-+\n",
        end_pattern=r"^\s*\-+\n",
        description="The title card of the Gaussian calculation task. link 101",
    )
    COORDS = MolOPPattern(
        end_pattern=r"(Basis read|Standard basis|Rotational constants \(GHZ\)|Symmetry turned off|The archive entry for this job was punched.)",
        content_pattern=r"\s+\d+\s+(\d+)\s+\d+\s+(\s*-?\d+\.\d+){3}",
        content_repeat=-1,
        description="The coordinates of the atoms in specific frame. link 202",
    )
    ROTATIONAL_CONST = MolOPPattern(
        content_pattern=r"Rotational constants \(GHZ\):\s*(\d+.\d+|\*+)\s*(\d+.\d+|\*+)\s*(\d+.\d+|\*+)",
        description="The rotational constants of the coordinates. link 202",
    )
    BASIS_SET = MolOPPattern(
        content_pattern=r"Standard basis:\s+([^\s]+)",
        description="The basis set used in the Gaussian calculation. link 301",
    )
    PSEUDOPOTENTIAL_PARAMETERS = MolOPPattern(
        start_pattern=r"^\s*={2,}\n",
        end_pattern=r"^\s*={2,}\n",
        end_offset=2,
        description="The parameters of the pseudopotential used in the Gaussian calculation. link 301",
    )
    DISPERSION_CORRECTION = MolOPPattern(
        content_pattern=r"^\s*R6Disp:\s*(.+)\s*Dispersion energy=(\s*-?\d+\.\d+)\s*Hartrees.",
        description="The dispersion correction used in the Gaussian calculation. link 301",
    )
    SOLVENT_PARAMETERS = MolOPPattern(
        start_pattern=r"Polarizable Continuum Model \(PCM\)",
        end_pattern=r"^\s*-{2,}\n",
        end_offset=2,
        description="The solvent parameters used in the Gaussian calculation. link 301",
    )
    SOLVENT_MODEL = MolOPPattern(
        content_pattern=r"Model\s*:\s*(.+)\.",
        description="The solvent model used in the Gaussian calculation. link 301",
    )
    SOLVENT_ATOM_RADII = MolOPPattern(
        content_pattern=r"Atomic radii\s*:\s*(.+)\.",
        description="The radii of the atoms in the solvent model. link 301",
    )
    SOLVENT_TYPE = MolOPPattern(
        content_pattern=r"Solvent\s*:\s*(.+),\s",
        description="The type of the solvent used in the Gaussian calculation. link 301",
    )
    SOLVENT_EPS = MolOPPattern(
        content_pattern=r"Eps\s*=\s*(\d+.\d*)",
        description="The solvent dielectric constant used in the Gaussian calculation. link 301",
    )
    SOLVENT_EPS_INF = MolOPPattern(
        content_pattern=r"Eps\((inf|infinity)\)\s*=\s*(\d+.\d*)",
        description="The solvent dielectric constant at infinity used in the Gaussian calculation. link 301",
    )
    SCF_ENERGY_AND_FUNCTIONAL = MolOPPattern(
        content_pattern=r"SCF Done:\s*E\(.*\)\s*=\s*(\s*-?\d+\.\d+)",
        description="The SCF energy and functional used in the Gaussian calculation. link 502",
    )
    SPIN_SPIN_SQUERE = MolOPPattern(
        content_pattern=r"<S\*\*2>=(\s*-?\d+\.\d+)\s+S=(\s*-?\d+\.\d+)",
        description="The total spin and spin squere exactly after the SCF calculation. link 502",
    )
    ENERGY_MP2_4 = MolOPPattern(
        content_pattern=r"E\d[\s\(\)SDTQ]*=\s*-*\d+.\d+D[+-]\d+\s*E*U(MP\d)[\(\)SDTQ]*\s*=\s*(-?\d+.\d+D[+-]\d+)",
        content_repeat=-1,
        description="The MP2-4 energy of the Gaussian calculation. link 804 for mp2 and link 913 for mp3-4",
    )
    ENERGY_MP5 = MolOPPattern(
        content_pattern=r"MP5\s*=\s*(-*\d+.\d+D[+-]\d+)\s*MP5\s*=\s*(-?\d+.\d+D[+-]\d+)",
        description="The MP5 energy of the Gaussian calculation. link 913",
    )
    ENERGY_CCSD = MolOPPattern(
        content_pattern=r"Wavefunction amplitudes converged. E\(Corr\)=\s*(-*\d+.\d*)",
        description="The CCSD energy of the Gaussian calculation. link 913",
    )
    ENERGY_CCSD_T = MolOPPattern(
        content_pattern=r"CCSD\(T\)\s*=\s*(-?\d+.\d+D[+-]\d+)",
        description="The CCSD(T) energy of the Gaussian calculation. link 913",
    )
    ISOTROPIC_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Isotropic polarizability for W=\s*\d+.\d+\s*(\d+.\d+)\s*Bohr\*\*3",
        description="The isotropic polarizability of the Gaussian calculation. link 1002",
    )
    POPULATION_ANALYSIS = MolOPPattern(
        start_pattern=r"Population analysis using the SCF [Dd]ensity.",
        end_pattern=r"^\s*N-N=.*",
        description="The population analysis of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS = MolOPPattern(
        content_pattern=r"(Alpha|Beta)\s*(occ.|virt.)\s*eigenvalues -- (.*)",
        content_repeat=-1,
        description="The molecular orbitals of the Gaussian calculation. link 601",
    )
    MULLIKEN_POPULATION = MolOPPattern(
        start_pattern=r"(Mulliken charges:|Mulliken atomic charges|Mulliken charges and spin densities:)",
        end_pattern=r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)",
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d+)",
        content_repeat=-1,
        description="The Mulliken population analysis of the Gaussian calculation. link 601",
    )
    MULLIKEN_SPIN_DENSITY = MolOPPattern(
        start_pattern=r"Mulliken charges and spin densities:",
        end_pattern=r"(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)",
        content_pattern=r"\d+\s+[A-Z][a-z]?(\s*-?\d+\.\d+){2}",
        content_repeat=-1,
        description="The Mulliken spin density analysis of the Gaussian calculation. link 601",
    )
    APT_POPULATION = MolOPPattern(
        start_pattern=r"APT\s*charges:",
        end_pattern=r"Sum of APT charges",
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d+)",
        content_repeat=-1,
        description="The APT population analysis of the Gaussian calculation. link 601",
    )
    LOWDIN_POPULATION = MolOPPattern(
        start_pattern=r"Lowdin\s*charges",
        end_pattern=r"Sum of Lowdin charges",
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d+)",
        content_repeat=-1,
        description="The Lowdin population analysis of the Gaussian calculation. link 601",
    )
    ELECTRONIC_SPATIAL_EXTENT = MolOPPattern(
        content_pattern=r"Electronic spatial extent \(au\):  <R\*\*2>=\s*(\d+\.\d+)",
        description="The electronic spatial extent (<R**2>) of the Gaussian calculation. link 601",
    )
    DIPOLE_MOMENT = MolOPPattern(
        start_pattern=r"Dipole moment \(field-independent basis",
        content_pattern=r"[XYZ]=\s+(\s*-?\d+\.\d+)",
        content_repeat=3,
        description="The dipole moment of the Gaussian calculation. link 601",
    )
    QUADRUPOLE_MOMENT = MolOPPattern(
        start_pattern=r"Quadrupole moment \(field-independent basis",
        content_pattern=r"\s[XYZ]{2}=\s+(\s*-?\d+\.\d+)",
        content_repeat=6,
        description="The quadrupole moment of the Gaussian calculation. link 601",
    )
    TRACELESS_QUADRUPOLE_MOMENT = MolOPPattern(
        start_pattern=r"Traceless Quadrupole moment \(field-independent basis",
        content_pattern=r"\s[XYZ]{2}=\s+(\s*-?\d+\.\d+)",
        content_repeat=6,
        description="The traceless quadrupole moment of the Gaussian calculation. link 601",
    )
    OCTAPOLE_MOMENT = MolOPPattern(
        start_pattern=r"Octapole moment \(field-independent basis",
        content_pattern=r"\s[XYZ]{3}=\s+(\s*-?\d+\.\d+)",
        content_repeat=10,
        description="The octapole moment of the Gaussian calculation. link 601",
    )
    HEXADECAPOLE_MOMENT = MolOPPattern(
        start_pattern=r"Hexadecapole moment \(field-independent basis",
        content_pattern=r"\s[XYZ]{4}=\s+(\s*-?\d+\.\d+)",
        content_repeat=15,
        description="The hexadecapole moment of the Gaussian calculation. link 601",
    )
    HIRSHFELD_POPULATION = MolOPPattern(
        start_pattern=r"Hirshfeld charges, spin densities, dipoles, and CM5 charges",
        end_pattern=r"Hirshfeld charges( and spin densities|) with hydrogens",
        content_pattern=r"^\s*\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d+){6}",
        content_repeat=-1,
        description="The Hirshfeld population analysis of the Gaussian calculation. link 601",
    )
    EXACT_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Exact\s*polarizability:(\s*-?\d+\.\d+){6}",
        description="The exact polarizability of the Gaussian calculation. link 601",
    )
    APPROX_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Approx\s*polarizability:(\s*-?\d+\.\d+){6}",
        description="The approximate polarizability of the Gaussian calculation. link 601",
    )
    ISOTROPIC_FERMI_CONTACT_COUPLING = MolOPPattern(
        start_pattern=r"Isotropic Fermi Contact Couplings",
        content_pattern=r"\s*\d+\s*[A-Z][a-z]?\((\d+)\)(\s*-?\d+\.\d+){4}",
        description="The isotropic Fermi contact coupling of the Gaussian calculation. link 601",
    )
    ESP_POPULATION = MolOPPattern(
        start_pattern=r"ESP charges:",
        end_pattern=r"Sum of ESP charges",
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(\s*-?\d+\.\d+)",
        content_repeat=-1,
        description="The ESP population analysis of the Gaussian calculation. link 601",
    )
    # TODO: add Bond order analysis patterns, link 607
    # TODO: add TDDFT patterns, link 914
    TDDFT_ORBITALS = MolOPPattern(description="TODO, link 914")
    FREQUENCY_ANALYSIS = MolOPPattern(
        end_pattern=r"^\s*(Grad){2,}",
        description="The frequency analysis of the Gaussian calculation. link 716",
    )
    DIPOLE_IN_FREQUENCY_ANALYSIS = MolOPPattern(
        content_pattern=r"^\s*Dipole\s*=(\s*-?\d*\.\d*D[+-]\d*){3}",
        description="The dipole showed again in the frequency analysis of the Gaussian calculation. link 716",
    )
    POLARIZIABILITIES_IN_FREQUENCY_ANALYSIS = MolOPPattern(
        content_pattern=r"^\s*Polarizability\s*=(\s*-?\d*\.\d*D[+-]\d*\n*){6}",
        description="The polarizabilities showed again in the frequency analysis of the Gaussian calculation. link 716",
    )
    DIAGONAL_VIBRATIONAL_POLARIZABILITY = MolOPPattern(
        content_pattern=r"^\s*Diagonal vibrational polarizability:\n(\s*-?\d*\.\d*){3}",
        description="The diagonal vibrational polarizability of the Gaussian calculation. link 716",
    )
    FREQUENCIES = MolOPPattern(
        content_pattern=r"Frequencies -- (\s+-?\d+\.\d+){1,3}",
        content_repeat=-1,
        description="The frequencies of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_REDUCED_MASS = MolOPPattern(
        content_pattern=r"Red. masses -- (\s+-?\d+\.\d+){1,3}",
        content_repeat=-1,
        description="The reduced masses of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_FORCE_CONSTANTS = MolOPPattern(
        content_pattern=r"Frc consts  -- (\s+-?\d+\.\d+){1,3}",
        content_repeat=-1,
        description="The force constants of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_IR_INTENSITIES = MolOPPattern(
        content_pattern=r"IR Inten    -- (\s+-?\d+\.\d+){1,3}",
        content_repeat=-1,
        description="The IR intensities of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_MODE = MolOPPattern(
        end_pattern=r"Thermochemistry",
        content_pattern=r"^\s+\d+\s+\d+\s*(\s*-?\d+\.\d+){3,9}",
        content_repeat=-1,
        description="The mode of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    TEMPEREATURE_PRESSURE = MolOPPattern(
        content_pattern=r"^\s*Temperature\s*(\d+\.\d+)\s*Kelvin\.\s*Pressure\s*(\d+\.\d+)\s*Atm\.",
        description="The temperature and pressure of the Gaussian calculation. link 716",
    )
    ROTATIONAL_TEMPERATURE = MolOPPattern(
        content_pattern=r"^\s*Rotational temperatures \(Kelvin\)(\s*\d+\.\d+){3}",
        description="The rotational temperatures of the Gaussian calculation. link 716",
    )
    ROTATIONAL_CONST_IN_FREQUENCY_ANALYSIS = MolOPPattern(
        content_pattern=r"^\s*Rotational constants \(GHZ\):(\s*\d+\.\d+){3}",
        description="The rotational constants of the Gaussian calculation. link 716",
    )
    VIBRATIONAL_TEMPERATURE = MolOPPattern(
        content_pattern=r"^\s*(Vibrational temperatures:|\(Kelvin\)|)(\s*\d+\.\d+){1,5}\n",
        content_repeat=-1,
        description="The vibrational temperatures of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_CORRECTION = MolOPPattern(
        content_pattern=r"(Zero-point|Thermal) correction(.*)=\s*(-?\d+\.\d*)",
        description="The thermochemistry corrections of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_SUM = MolOPPattern(
        content_pattern=r"Sum of electronic and (thermal Free Energies|thermal Enthalpies|thermal Energies|zero-point Energies)=\s*(-?\d+.\d*)",
        description="The thermochemistry sums of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_CV_S = MolOPPattern(
        start_pattern=r"E \(Thermal\)\s*CV\s*S",
        content_pattern=r"^\sTotal(\s*-?\d+\.\d+){3}",
        description="The thermochemistry CV and S of the Gaussian calculation. link 716",
    )
    FORCES_IN_CARTESIAN = MolOPPattern(
        start_pattern=r"^\s*Center\s+Atomic\s+Forces \(Hartrees/Bohr\)\s*",
        end_pattern=r"^\s*(-){2,}",
        end_offset=-1,
        content_pattern=r"^\s*\d+\s+\d+(\s*-?\d+\.\d+){3}",
        description="The forces in Cartesian coordinates of the Gaussian calculation. link 716",
    )
    HESSIAN_IN_CARTESIAN = MolOPPattern(
        start_pattern=r"^\s*Force constants in Cartesian coordinates",
        end_pattern=r"^\s*FormGI is forming",
        end_offset=-1,
        content_pattern=r"-?\d+.\d+D[+-]\d+",
        description="The Hessian in Cartesian coordinates of the Gaussian calculation. link 716",
    )
    BERNY_STATE_MAJOR = MolOPPattern(
        start_pattern=r"^\s*\(Enter\s*.*103.exe\)\n",
        end_pattern=r"^\s*Leave Link  103.*MaxMem=\s*(\d+)\s*cpu:\s*(\d+\.\d+)",
        content_pattern=r"(Maximum\s+Force|RMS\s+Force|Maximum\s+Displacement|RMS\s+Displacement)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(NO|YES)",
        content_repeat=4,
        description="The Berny optimization state. Patterns showed when #p used. link 103",
    )
    BERNY_STATE_BACKUP = MolOPPattern(
        start_pattern=r"^\s*(Grad){2,}",
        end_pattern=r"^\s*(Grad){2,}",
        content_pattern=r"(Maximum\s+Force|RMS\s+Force|Maximum\s+Displacement|RMS\s+Displacement)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(NO|YES)",
        content_repeat=4,
        description="The Berny optimization state. link 103",
    )
    ELECTRIC_DIPOLE_MOMENT = MolOPPattern(
        start_pattern=r"^\s*Electric dipole moment \(input orientation\):",
        content_pattern=r"^\s*(Tot|x|y|z)(\s*-?\d+\.\d+D[+-]\d+){3}",
        content_repeat=4,
        description="The electric dipole moment of the Gaussian calculation. link 9999",
    )
    DIPOLE_POLARIZABILITY = MolOPPattern(
        start_pattern=r"^\s*Dipole polarizability, Alpha \(input orientation\)",
        content_pattern=r"^\s*(iso|aniso|xx|yx|yy|zx|zy|zz)(\s*-?\d+\.\d+D[+-]\d+){3}",
        content_repeat=8,
        description="The dipole polarizability of the Gaussian calculation. link 9999",
    )
    ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"^\s*(1\\1\\GINC|1\|1\|GINC)",
        end_pattern=r"(\\@\n|@\n)",
        description="The archive tail of the Gaussian calculation. link 9999",
    )
    JOB_TYPE_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        content_pattern=r".*\\",
        description="The job type in the 4th block of the archive tail of the Gaussian calculation. link 9999",
    )
    FUNCTIONAL_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        content_pattern=r".*\\",
        description="The functional in the 5th block of the archive tail of the Gaussian calculation. link 9999",
    )
    BASIS_SET_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        content_pattern=r".*\\",
        description="The basis set in the 6th block of the archive tail of the Gaussian calculation. link 9999",
    )
    KEYWORDS_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        start_offset=5,
        content_pattern=r".*\\",
        description="The keywords in the 12th block of the archive tail of the Gaussian calculation. link 9999",
    )
    TITLE_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        start_offset=1,
        content_pattern=r".*\\",
        description="The title in the 14th block of the archive tail of the Gaussian calculation. link 9999",
    )
    CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        start_offset=1,
        content_pattern=r".*\\",
        description="The charge and spin multiplicity in the 16th block of the archive tail of the Gaussian calculation. link 9999",
    )
    COORS_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"\\([A-Z][a-z]?),(-?\d+\.\d+),(-?\d+\.\d+),(-?\d+\.\d+)",
        content_repeat=-1,
        description="The coordinates in the 17th block of the archive tail of the Gaussian calculation. link 9999",
    )
    VERSION_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"\\",
        start_offset=1,
        content_pattern=r".*\\",
        description="The version in the archive tail of the Gaussian calculation. link 9999",
    )
    ENERGIES_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"(HF|MP2|MP3|MP4[SDTQ]*|CCSD[\(\)T]*)=(\s*-?\d+\.\d+)",
        content_repeat=-1,
        description="The energies in the archive tail of the Gaussian calculation. link 9999",
    )
    THERMOCHEMISTRY_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"(ZeroPoint|Thermal|ETot|HTot|GTot)=(\s*-?\d+\.\d+)",
        description="The thermochemistry in the archive tail of the Gaussian calculation. link 9999",
    )
    DIPOLE_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Dipole\s*=(\s*-?\d+\.\d+){3}",
        description="The dipole in the archive tail of the Gaussian calculation. link 9999",
    )
    POLAR_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Polar\s*=\s*(-?\d+\.\d+){6}",
        description="The polarizability in the archive tail of the Gaussian calculation. link 9999",
    )
    QUADRUPOLE_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Quadrupole\s*=\s*(-?\d+\.\d+){6}",
        description="The quadrupole in the archive tail of the Gaussian calculation. link 9999",
    )
    HESSIAN_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"NImag=\d+\\\\",
        content_pattern=r"(-?\d+\.\d+)",
        content_repeat=-1,
        description="The Hessian in the archive tail of the Gaussian calculation. link 9999",
    )
    JOB_CPU_TIME = MolOPPattern(
        content_pattern=r"^\s*Job cpu time:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes\s*(\d+\.\d+)\s*seconds",
        description="The job cpu time of the Gaussian calculation. link 9999",
    )
    ELAPSED_TIME = MolOPPattern(
        content_pattern=r"^\s*Elapsed time:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes\s*(\d+\.\d+)\s*seconds",
        description="The elapsed time of the Gaussian calculation. link 9999",
    )
    TERMINATION_STATUS = MolOPPattern(
        content_pattern=r"^\s*(Normal|Error) termination (.+)\s*(Mon|Tue|Wed|Thu|Fri|Sat|Sun)\s*(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s*(\d+)\s*(\d+)\s*:\s*(\d+)\s*:\s*(\d+)\s*(\d+)",
        description="The termination status of the Gaussian calculation. link 9999",
    )


class G16FchkPatterns:
    """
    G16 fchk file patterns
    """
