"""
Author: TMJ
Date: 2025-02-16 19:20:04
LastEditors: TMJ
LastEditTime: 2026-03-23 15:56:46
Description: 请填写简介
"""

from molop.io.base_models.SearchPattern import MolOPPattern


_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][-+]?\d+)?"
_ROTATIONAL_CONSTANT = rf"(?:{_FLOAT}|\*+)"


class G16LogPatterns:
    """
    G16 log file patterns
    """

    FLOAT_TOKEN = MolOPPattern(
        content_pattern=r"(?P<value>[-+]?\d+\.\d+(?:[DEde][-+]?\d+)?)",
        content_repeat=0,
    )
    FLOAT_TOKEN_5DP = MolOPPattern(
        content_pattern=r"(?P<value>[-+]?\d+\.\d{5}(?:[DEde][-+]?\d+)?)",
        content_repeat=0,
    )

    def float_token(self, decimal_places: int | None = None) -> MolOPPattern:
        if decimal_places is None:
            return self.FLOAT_TOKEN
        if decimal_places == 5:
            return self.FLOAT_TOKEN_5DP
        return MolOPPattern(
            content_pattern=rf"(?P<value>[-+]?\d+\.\d{{{decimal_places}}}(?:[DEde][-+]?\d+)?)",
            content_repeat=0,
        )

    LINK1_SECTION = MolOPPattern(
        content_pattern=r"^\s*(?:Link1:\s+Proceeding to internal job step number\s+\d+\.|"
        r"Entering Link 1 = .*?)\s*$",
        content_repeat=0,
    )
    VERSION = MolOPPattern(
        start_pattern=r"^\s\*+\n",
        end_pattern=r"^\s\*+\n",
        end_offset=1,
        content_pattern=r"(?P<version>Gaussian\s+\d+:\s*[^\s]*\s*"
        r"\d+-[A-Za-z]{3}-\d{4}\n\s+\d+-[A-Za-z]{3}-\d{4})",
        description="The exact version of Gaussian used. link 1",
    )
    VERSION_TOKEN = MolOPPattern(
        content_pattern=r"\b(?P<version>[A-Za-z0-9]+-G16Rev[A-Za-z0-9.]+)\b",
        description="The normalized Gaussian 16 version token from the banner.",
    )
    OPTIONS = MolOPPattern(
        start_pattern=r"^\s\*+\n",
        end_pattern=r"^\s*-+\n",
        content_pattern=r"^\s*(?P<key>%[a-zA-Z0-9]+)=(?P<value>[^\s]+)",
        content_repeat=0,
        description="The options used in the Gaussian calculation. Such as: MEMORY, CPU, etc. link 1",
    )
    KEYWORDS = MolOPPattern(
        start_pattern=r"^\s*-+\n",
        end_pattern=r"^\s*-+\n",
        end_offset=1,
        content_pattern=r"^\s*(?P<keywords>#[a-zA-Z0-9\(\)=,\/\\\s\n\-]+)\n\s*\-+",
        description="The keywords used in the Gaussian calculation. link 1",
    )
    TITLE = MolOPPattern(
        start_pattern=r"^\s*-+\n",
        start_offset=1,
        end_pattern=r"^\s*-+\n",
        end_offset=1,
        content_pattern=r"----\n\s*(?P<title>(?:.\n*)+)\n\s*----",
        description="The title card of the Gaussian calculation task. link 101",
    )
    CHARGE_MULTIPLICITY = MolOPPattern(
        content_pattern=r"Charge\s*=\s*(?P<charge>[-\+\d]+)\s+"
        r"Multiplicity\s*=\s*(?P<multiplicity>\d+)",
        description="The charge and multiplicity of the Gaussian calculation. link 101",
    )
    INITIAL_INPUT_COORDS = MolOPPattern(
        start_pattern="Symbolic Z-matrix:",
        start_regex=False,
        end_pattern=r"^\s*\n",
        content_pattern=r"^\s*(?P<symbol>[A-Za-z]{1,3})\s+(?P<x>\s*-?\d+\.\d*)"
        r"(?P<y>\s*-?\d+\.\d*)(?P<z>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The initial input coordinates of the Gaussian calculation. link 101",
    )
    INPUT_COORDS = MolOPPattern(
        start_pattern="Input orientation:",
        start_regex=False,
        end_pattern=r"^\s*-+\n",
        end_offset=2,
        content_pattern=r"\s+\d+\s+(?P<atomic_number>\d+)\s+\d+\s+(?P<x>\s*-?\d+\.\d*)"
        r"(?P<y>\s*-?\d+\.\d*)(?P<z>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The input coordinates of the Gaussian calculation. link 202",
    )
    STANDARD_COORDS = MolOPPattern(
        start_pattern="Standard orientation:",
        start_regex=False,
        end_pattern=r"^\s*-+\n",
        end_offset=2,
        content_pattern=r"\s+\d+\s+(?P<atomic_number>\d+)\s+\d+\s+(?P<x>\s*-?\d+\.\d*)"
        r"(?P<y>\s*-?\d+\.\d*)(?P<z>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The standard coordinates of the Gaussian calculation. link 202",
    )
    COORDS = MolOPPattern(
        end_pattern=r"(?:Basis read|Standard basis|Rotational constants \(GHZ\)|Symmetry turned off|The archive entry for this job was punched.)",
        content_pattern=r"^\s*\d+\s+(?P<atomic_number>\d+)\s+\d+\s+(?P<x>\s*-?\d+\.\d*)"
        r"(?P<y>\s*-?\d+\.\d*)(?P<z>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The coordinates of the atoms in specific frame. link 202",
    )
    ROTATIONAL_CONST = MolOPPattern(
        content_pattern=r"^\s*Rotational constants \(GHZ\):\s*"
        rf"(?P<a>{_ROTATIONAL_CONSTANT})\s+"
        rf"(?P<b>{_ROTATIONAL_CONSTANT})\s+"
        rf"(?P<c>{_ROTATIONAL_CONSTANT})",
        description="The rotational constants of the coordinates. link 202",
    )
    BASIS_SET = MolOPPattern(
        content_pattern=r"^\s*Standard basis:\s+(?P<basis>[^\s]+)",
        description="The basis set used in the Gaussian calculation. link 301",
    )
    PSEUDOPOTENTIAL_PARAMETERS = MolOPPattern(
        start_pattern=r"^\s*={2,}\n",
        end_pattern=r"^\s*={2,}\n",
        end_offset=2,
        description="The parameters of the pseudopotential used in the Gaussian calculation. link 301",
    )
    DISPERSION_CORRECTION = MolOPPattern(
        content_pattern=r"^\s*R6Disp:\s*(?P<label>.+)\s*Dispersion energy="
        r"(?P<energy>\s*-?\d+\.\d*)\s*Hartrees.",
        description="The dispersion correction used in the Gaussian calculation. link 301",
    )
    SOLVENT_PARAMETERS = MolOPPattern(
        start_pattern="Polarizable Continuum Model (PCM)",
        start_regex=False,
        end_pattern=r"^\s*-{2,}\n",
        description="The solvent parameters used in the Gaussian calculation. link 301",
    )
    SOLVENT_MODEL = MolOPPattern(
        content_pattern=r"Model\s*:\s*(?P<model>.+)\.",
        description="The solvent model used in the Gaussian calculation. link 301",
    )
    SOLVENT_ATOM_RADII = MolOPPattern(
        content_pattern=r"Atomic radii\s*:\s*(?P<radii>.+)\.",
        description="The radii of the atoms in the solvent model. link 301",
    )
    SOLVENT_TYPE = MolOPPattern(
        content_pattern=r"Solvent\s*:\s*(?P<solvent>.+),\s",
        description="The type of the solvent used in the Gaussian calculation. link 301",
    )
    SOLVENT_EPS = MolOPPattern(
        content_pattern=r"Eps\s*=\s*(?P<value>\d+.\d*)",
        description="The solvent dielectric constant used in the Gaussian calculation. link 301",
    )
    SOLVENT_EPS_INF = MolOPPattern(
        content_pattern=r"Eps\((?P<label>inf|infinity)\)\s*=\s*(?P<value>\d+.\d*)",
        description="The solvent dielectric constant at infinity used in the Gaussian calculation. link 301",
    )
    SCF_ENERGIES = MolOPPattern(
        start_pattern=r"^\s*SCF Done:",
        description="The SCF energies of the Gaussian calculation. link 502",
    )
    SCF_ENERGY_AND_FUNCTIONAL = MolOPPattern(
        content_pattern=r"^\s*SCF Done:\s*E\((?P<functional>.*)\)\s*=\s*"
        r"(?P<energy>\s*-?\d+\.\d*)",
        description="The SCF energy and functional used in the Gaussian calculation. link 502",
    )
    EXTERNAL_ENERGY_RESULT = MolOPPattern(
        content_pattern=(
            r"^[ \t]*External calculation of energy"
            r"(?: and first derivatives|, first and second derivatives)?\.[ \t]*\r?\n"
            r"(?:(?![ \t]*(?:External calculation of energy|Error termination|"
            r"Normal termination))[^\r\n]*(?:\r?\n|$)){0,32}?"
            rf"^[ \t]*Energy=[ \t]*(?P<energy>{_FLOAT})[ \t]+"
            r"NIter=[ \t]*(?P<niter>\d+)\.[ \t]*$"
        ),
        content_repeat=0,
        description=(
            "An energy returned successfully through Gaussian's External interface. link 402"
        ),
    )
    SPIN_SPIN_SQUERE = MolOPPattern(
        content_pattern=r"<S\*\*2>=(?P<spin_square>\s*-?\d+\.\d*)\s+"
        r"S=(?P<spin_quantum_number>\s*-?\d+\.\d*)",
        description="The total spin and spin squere exactly after the SCF calculation. link 502",
    )
    ENERGY_MP2_4 = MolOPPattern(
        content_pattern=r"E\d[\s\(\)SDTQ]*=\s*-*\d+.\d+D[+-]\d+\s*"
        r"E*U(?P<method>MP\d)[\(\)SDTQ]*\s*=\s*(?P<energy>-?\d+.\d+D[+-]\d+)",
        content_repeat=0,
        description="The MP2-4 energy of the Gaussian calculation. link 804 for mp2 and link 913 for mp3-4",
    )
    ENERGY_MP5 = MolOPPattern(
        content_pattern=r"MP5\s*=\s*(?P<first_energy>-*\d+.\d+D[+-]\d+)\s*"
        r"MP5\s*=\s*(?P<energy>-?\d+.\d+D[+-]\d+)",
        description="The MP5 energy of the Gaussian calculation. link 913",
    )
    ENERGY_CCSD = MolOPPattern(
        content_pattern=r"Wavefunction amplitudes converged. E\(Corr\)=\s*(?P<energy>-*\d+.\d*)",
        description="The CCSD energy of the Gaussian calculation. link 913",
    )
    ENERGY_CCSD_T = MolOPPattern(
        content_pattern=r"CCSD\(T\)\s*=\s*(?P<energy>-?\d+.\d+D[+-]\d+)",
        description="The CCSD(T) energy of the Gaussian calculation. link 913",
    )
    ISOTROPIC_POLARIZABILITY = MolOPPattern(
        start_pattern="Isotropic polarizability for W=",
        start_regex=False,
        end_pattern=r"^\s*Isotropic polarizability for W=\s*\d+.\d+\s*(?:\d+.\d+)\s*Bohr\*\*3",
        content_pattern=r"\d+.\d+\s*(?P<value>\d+.\d+)",
        description="The isotropic polarizability of the Gaussian calculation. link 1002",
    )
    NMR_SHIELDING = MolOPPattern(
        start_pattern=r"^\s*[A-Za-z0-9()\-+]+\s+[A-Za-z0-9()\-+]+\s+"
        r"Magnetic shielding tensor \(ppm\):\s*$",
        end_pattern=r"^\s*End of Minotr F\.D\. properties file",
        content_pattern=(
            rf"^\s*(?P<atom_index>\d+)\s+(?P<atom_symbol>[A-Z][a-z]?)\s+"
            rf"Isotropic\s*=\s*(?P<isotropic>{_FLOAT})\s+"
            rf"Anisotropy\s*=\s*(?P<anisotropy>{_FLOAT})\s*$\n"
            rf"\s*XX=\s*(?P<xx>{_FLOAT})\s+YX=\s*(?P<yx>{_FLOAT})\s+"
            rf"ZX=\s*(?P<zx>{_FLOAT})\s*$\n"
            rf"\s*XY=\s*(?P<xy>{_FLOAT})\s+YY=\s*(?P<yy>{_FLOAT})\s+"
            rf"ZY=\s*(?P<zy>{_FLOAT})\s*$\n"
            rf"\s*XZ=\s*(?P<xz>{_FLOAT})\s+YZ=\s*(?P<yz>{_FLOAT})\s+"
            rf"ZZ=\s*(?P<zz>{_FLOAT})\s*$\n"
            rf"\s*Eigenvalues:\s*(?P<eigenvalue_1>{_FLOAT})\s+"
            rf"(?P<eigenvalue_2>{_FLOAT})\s+(?P<eigenvalue_3>{_FLOAT})\s*$"
        ),
        content_repeat=0,
        description="Per-atom Gaussian magnetic shielding tensors. link 1002",
    )
    NMR_SHIELDING_HEADER = MolOPPattern(
        content_pattern=r"^\s*(?P<method>[A-Za-z0-9()\-+]+)\s+"
        r"(?P<gauge>[A-Za-z0-9()\-+]+)\s+Magnetic shielding tensor \(ppm\):\s*$",
        description="Gaussian NMR shielding method and gauge header. link 1002",
    )
    NMR_COUPLING_FC_K = MolOPPattern(
        start_pattern=r"^\s*Fermi Contact \(FC\) contribution to K \(Hz\):\s*$",
        end_pattern=r"^\s*Fermi Contact \(FC\) contribution to J \(Hz\):\s*$",
        description="Gaussian Fermi-contact reduced spin-spin coupling matrix. link 1002",
    )
    NMR_COUPLING_FC_J = MolOPPattern(
        start_pattern=r"^\s*Fermi Contact \(FC\) contribution to J \(Hz\):\s*$",
        end_pattern=r"^\s*Spin-dipolar \(SD\) contribution to K \(Hz\):\s*$",
        description="Gaussian Fermi-contact spin-spin coupling matrix. link 1002",
    )
    NMR_COUPLING_SD_K = MolOPPattern(
        start_pattern=r"^\s*Spin-dipolar \(SD\) contribution to K \(Hz\):\s*$",
        end_pattern=r"^\s*Spin-dipolar \(SD\) contribution to J \(Hz\):\s*$",
        description="Gaussian spin-dipolar reduced spin-spin coupling matrix. link 1002",
    )
    NMR_COUPLING_SD_J = MolOPPattern(
        start_pattern=r"^\s*Spin-dipolar \(SD\) contribution to J \(Hz\):\s*$",
        end_pattern=r"^\s*Paramagnetic spin-orbit \(PSO\) contribution to K \(Hz\):\s*$",
        description="Gaussian spin-dipolar spin-spin coupling matrix. link 1002",
    )
    NMR_COUPLING_PSO_K = MolOPPattern(
        start_pattern=r"^\s*Paramagnetic spin-orbit \(PSO\) contribution to K \(Hz\):\s*$",
        end_pattern=r"^\s*Paramagnetic spin-orbit \(PSO\) contribution to J \(Hz\):\s*$",
        description="Gaussian paramagnetic spin-orbit reduced coupling matrix. link 1002",
    )
    NMR_COUPLING_PSO_J = MolOPPattern(
        start_pattern=r"^\s*Paramagnetic spin-orbit \(PSO\) contribution to J \(Hz\):\s*$",
        end_pattern=r"^\s*Diamagnetic spin-orbit \(DSO\) contribution to K \(Hz\):\s*$",
        description="Gaussian paramagnetic spin-orbit coupling matrix. link 1002",
    )
    NMR_COUPLING_DSO_K = MolOPPattern(
        start_pattern=r"^\s*Diamagnetic spin-orbit \(DSO\) contribution to K \(Hz\):\s*$",
        end_pattern=r"^\s*Diamagnetic spin-orbit \(DSO\) contribution to J \(Hz\):\s*$",
        description="Gaussian diamagnetic spin-orbit reduced coupling matrix. link 1002",
    )
    NMR_COUPLING_DSO_J = MolOPPattern(
        start_pattern=r"^\s*Diamagnetic spin-orbit \(DSO\) contribution to J \(Hz\):\s*$",
        end_pattern=r"^\s*Total nuclear spin-spin coupling K \(Hz\):\s*$",
        description="Gaussian diamagnetic spin-orbit coupling matrix. link 1002",
    )
    NMR_TOTAL_COUPLING_K = MolOPPattern(
        start_pattern=r"^\s*Total nuclear spin-spin coupling K \(Hz\):\s*$",
        end_pattern=r"^\s*Total nuclear spin-spin coupling J \(Hz\):\s*$",
        description="Gaussian total reduced spin-spin coupling matrix. link 1002",
    )
    NMR_TOTAL_COUPLING_J = MolOPPattern(
        start_pattern=r"^\s*Total nuclear spin-spin coupling J \(Hz\):\s*$",
        end_pattern=r"^\s*End of Minotr F\.D\. properties file",
        description="Gaussian total spin-spin coupling matrix. link 1002",
    )
    NMR_COUPLING_COLUMN_HEADER = MolOPPattern(
        content_pattern=r"^\s*(?P<columns>\d+(?:\s+\d+)*)\s*$",
        description="Gaussian spin-spin coupling matrix column labels.",
    )
    NMR_COUPLING_ROW = MolOPPattern(
        content_pattern=r"^\s*(?P<row>\d+)\s+"
        r"(?P<values>[-+]?\d+\.\d+(?:[DEde][-+]?\d+)"
        r"(?:\s+[-+]?\d+\.\d+(?:[DEde][-+]?\d+))*)\s*$",
        description="Gaussian spin-spin coupling matrix row.",
    )
    POPULATION_ANALYSIS = MolOPPattern(
        start_pattern=r"^\s*Population analysis using the (?:SCF|CC) [Dd]ensity.",
        end_pattern=r"^\s*N-N=.*",
        description="The population analysis of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS_SYMMETRY_ALPHA = MolOPPattern(
        start_pattern="Alpha Orbitals:",
        start_regex=False,
        end_pattern="Beta  Orbitals:",
        end_regex=False,
        content_pattern=r"^\s*(?P<occupancy>Occupied|Virtual|)\s*"
        r"(?P<symbols>(?:\(\??[A-Z0-9]+\)\s*){1,12})\n",
        content_repeat=0,
        description="The alpha molecular orbitals symmetry of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS_SYMMETRY_BETA = MolOPPattern(
        start_pattern="Beta  Orbitals:",
        start_regex=False,
        end_pattern=r"^\s*(?:The electronic state is (?:.*)\.|"
        r"(?:Alpha|Beta)\s+(?:occ\.|virt\.) eigenvalues --)",
        content_pattern=r"^\s*(?P<occupancy>Occupied|Virtual|)\s*"
        r"(?P<symbols>(?:\(\??[A-Z0-9]+\)\s*){1,12})\n",
        content_repeat=0,
        description="The beta molecular orbitals symmetry of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS_SYMMETRY = MolOPPattern(
        start_pattern="Orbital symmetries:",
        start_regex=False,
        end_pattern=r"^\s*The electronic state is (?:.*)\.",
        content_pattern=r"^\s*(?P<occupancy>Occupied|Virtual|)\s*"
        r"(?P<symbols>(?:\(\??[A-Z0-9]+\)\s*){1,12})\n",
        content_repeat=0,
        description="The molecular orbitals symmetry of the Gaussian calculation. link 601",
    )
    ELECTRONIC_STATE = MolOPPattern(
        content_pattern=r"^\s*The electronic state is (?P<state>.*)\.",
        description="The electronic state of the Gaussian calculation. link 601",
    )
    MOLECULAR_ORBITALS = MolOPPattern(
        content_pattern=r"^\s*(?P<spin>Alpha|Beta)\s*(?P<occupancy>occ.|virt.)\s*"
        r"eigenvalues -- (?P<values>.*)",
        content_repeat=0,
        description="The molecular orbitals of the Gaussian calculation. link 601",
    )
    MULLIKEN_POPULATION = MolOPPattern(
        start_pattern=r"(?:Mulliken charges:|Mulliken atomic charges|Mulliken charges and spin densities:)",
        end_pattern=r"(?:Sum of Mulliken )(?:.*)(?:charges)\s*=\s*(?:\D)",
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(?P<charge>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Mulliken population analysis of the Gaussian calculation. link 601",
    )
    MULLIKEN_SPIN_DENSITY = MolOPPattern(
        start_pattern="Mulliken charges and spin densities:",
        start_regex=False,
        end_pattern=r"(?:Sum of Mulliken )(?:.*)(?:charges)\s*=\s*(?:\D)",
        content_pattern=r"\d+\s+[A-Z][a-z]?(?P<charge>\s*-?\d+\.\d*)"
        r"(?P<spin>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Mulliken spin density analysis of the Gaussian calculation. link 601",
    )
    APT_POPULATION = MolOPPattern(
        start_pattern="APT charges:",
        start_regex=False,
        end_pattern="Sum of APT charges",
        end_regex=False,
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(?P<charge>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The APT population analysis of the Gaussian calculation. link 601",
    )
    LOWDIN_POPULATION = MolOPPattern(
        start_pattern="Lowdin charges",
        start_regex=False,
        end_pattern="Sum of Lowdin charges",
        end_regex=False,
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(?P<charge>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Lowdin population analysis of the Gaussian calculation. link 601",
    )
    ELECTRONIC_SPATIAL_EXTENT = MolOPPattern(
        content_pattern=r"Electronic spatial extent \(au\):  <R\*\*2>=\s*(?P<value>\d+\.\d*)",
        description="The electronic spatial extent (<R**2>) of the Gaussian calculation. link 601",
    )
    DIPOLE_MOMENT = MolOPPattern(
        start_pattern="Dipole moment (field-independent basis",
        start_regex=False,
        content_pattern=r"[XYZ]=\s+(?P<value>\s*-?\d+\.\d*)",
        content_repeat=3,
        description="The dipole moment of the Gaussian calculation. link 601",
    )
    QUADRUPOLE_MOMENT = MolOPPattern(
        start_pattern="Quadrupole moment (field-independent basis",
        start_regex=False,
        content_pattern=r"\s[XYZ]{2}=\s+(?P<value>\s*-?\d+\.\d*)",
        content_repeat=6,
        description="The quadrupole moment of the Gaussian calculation. link 601",
    )
    TRACELESS_QUADRUPOLE_MOMENT = MolOPPattern(
        start_pattern="Traceless Quadrupole moment (field-independent basis",
        start_regex=False,
        content_pattern=r"\s[XYZ]{2}=\s+(?P<value>\s*-?\d+\.\d*)",
        content_repeat=6,
        description="The traceless quadrupole moment of the Gaussian calculation. link 601",
    )
    OCTAPOLE_MOMENT = MolOPPattern(
        start_pattern="Octapole moment (field-independent basis",
        start_regex=False,
        content_pattern=r"\s[XYZ]{3}=\s+(?P<value>\s*-?\d+\.\d*)",
        content_repeat=10,
        description="The octapole moment of the Gaussian calculation. link 601",
    )
    HEXADECAPOLE_MOMENT = MolOPPattern(
        start_pattern="Hexadecapole moment (field-independent basis",
        start_regex=False,
        content_pattern=r"\s[XYZ]{4}=\s+(?P<value>\s*-?\d+\.\d*)",
        content_repeat=15,
        description="The hexadecapole moment of the Gaussian calculation. link 601",
    )
    HIRSHFELD_POPULATION = MolOPPattern(
        start_pattern="Hirshfeld charges, spin densities, dipoles, and CM5 charges",
        start_regex=False,
        end_pattern=r"^\s*Hirshfeld charges(?: and spin densities)? with hydrogens",
        content_pattern=r"^\s*\d+\s+[A-Z][a-z]?\s+(?P<charge>\s*-?\d+\.\d*)"
        r"(?P<spin>\s*-?\d+\.\d*)(?P<dipole_x>\s*-?\d+\.\d*)"
        r"(?P<dipole_y>\s*-?\d+\.\d*)(?P<dipole_z>\s*-?\d+\.\d*)"
        r"(?P<q_cm5>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The Hirshfeld population analysis of the Gaussian calculation. link 601",
    )
    EXACT_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Exact\s*polarizability:(?P<xx>\s*-?\d+\.\d*)"
        r"(?P<xy>\s*-?\d+\.\d*)(?P<yy>\s*-?\d+\.\d*)"
        r"(?P<xz>\s*-?\d+\.\d*)(?P<yz>\s*-?\d+\.\d*)"
        r"(?P<zz>\s*-?\d+\.\d*)",
        description="The exact polarizability of the Gaussian calculation. link 601",
    )
    APPROX_POLARIZABILITY = MolOPPattern(
        content_pattern=r"Approx\s*polarizability:(?P<xx>\s*-?\d+\.\d*)"
        r"(?P<xy>\s*-?\d+\.\d*)(?P<yy>\s*-?\d+\.\d*)"
        r"(?P<xz>\s*-?\d+\.\d*)(?P<yz>\s*-?\d+\.\d*)"
        r"(?P<zz>\s*-?\d+\.\d*)",
        description="The approximate polarizability of the Gaussian calculation. link 601",
    )
    ISOTROPIC_FERMI_CONTACT_COUPLING = MolOPPattern(
        start_pattern="Isotropic Fermi Contact Couplings",
        start_regex=False,
        content_pattern=r"\s*\d+\s*[A-Z][a-z]?\((?P<isotope>\d+)\)"
        r"(?P<isotropic>\s*-?\d+\.\d*)(?P<anisotropic>\s*-?\d+\.\d*)"
        r"(?P<dipolar>\s*-?\d+\.\d*)(?P<contact>\s*-?\d+\.\d*)",
        description="The isotropic Fermi contact coupling of the Gaussian calculation. link 601",
    )
    ESP_POPULATION = MolOPPattern(
        start_pattern="ESP charges:",
        start_regex=False,
        end_pattern="Sum of ESP charges",
        end_regex=False,
        content_pattern=r"\d+\s+[A-Z][a-z]?\s+(?P<charge>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The ESP population analysis of the Gaussian calculation. link 601",
    )
    NPA_POPULATION = MolOPPattern(
        start_pattern="Summary of Natural Population Analysis:",
        start_regex=False,
        end_pattern=r"^\s*=+\s*$",
        content_pattern=(
            r"^\s*[A-Z][a-z]?\s+\d+\s+"
            r"(?P<charge>[+-]?\d+\.\d+)\s+"
            r"[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+"
            r"[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s*$"
        ),
        content_repeat=0,
        description="The natural population analysis summary. NBO link 607",
    )
    DIPOLE_BEFORE_FORCE = MolOPPattern(
        start_pattern="Dipole        =",
        start_regex=False,
        content_pattern=r"(?P<value>\s*-*\d+.\d+D[+-]\d+)",
        content_repeat=3,
        description="The dipole moment before the force calculation. link 716",
    )
    POLARIZIABILITIES_BEFORE_FORCE = MolOPPattern(
        start_pattern="Polarizability=",
        start_regex=False,
        content_pattern=r"(?P<value>\s*-*\d+.\d+D[+-]\d+)",
        content_repeat=6,
        description="The polarizabilities before the force calculation. link 716",
    )
    # TODO: add Bond order analysis patterns, link 607
    # TODO: add TDDFT patterns, link 914
    TDDFT_ORBITALS = MolOPPattern(description="TODO, link 914")
    FREQUENCY_ANALYSIS = MolOPPattern(
        start_pattern="Harmonic frequencies (cm**-1)",
        start_regex=False,
        end_pattern=r" \n",
        end_regex=False,
        description="The frequency analysis of the Gaussian calculation. link 716",
    )
    DIAGONAL_VIBRATIONAL_POLARIZABILITY = MolOPPattern(
        content_pattern=r"^\s*Diagonal vibrational polarizability:\n"
        r"(?P<x>\s*-?\d*\.\d*)(?P<y>\s*-?\d*\.\d*)(?P<z>\s*-?\d*\.\d*)",
        description="The diagonal vibrational polarizability of the Gaussian calculation. link 716",
    )
    FREQUENCIES = MolOPPattern(
        content_pattern=r"Frequencies -- (?P<values>(?:\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The frequencies of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_REDUCED_MASS = MolOPPattern(
        content_pattern=r"Red. masses -- (?P<values>(?:\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The reduced masses of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_FORCE_CONSTANTS = MolOPPattern(
        content_pattern=r"Frc consts  -- (?P<values>(?:\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The force constants of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_IR_INTENSITIES = MolOPPattern(
        content_pattern=r"IR Inten    -- (?P<values>(?:\s+-?\d+\.\d*){1,3})",
        content_repeat=0,
        description="The IR intensities of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    FREQUENCIES_MODE = MolOPPattern(
        end_pattern="Thermochemistry",
        end_regex=False,
        content_pattern=r"^\s+\d+\s+\d+\s*(?P<values>(?:\s*-?\d+\.\d*){3,9})",
        content_repeat=0,
        description="The mode of each vibration in the frequency analysis of the Gaussian calculation. link 716",
    )
    TEMPEREATURE_PRESSURE = MolOPPattern(
        start_pattern="- Thermochemistry",
        start_regex=False,
        end_pattern=r"^\s*Temperature\s*\d+\.\d*\s*Kelvin\.\s*Pressure\s*\d+\.\d*\s*Atm\.",
        content_pattern=r"^\s*Temperature\s*(?P<temperature>\d+\.\d*)\s*Kelvin\.\s*"
        r"Pressure\s*(?P<pressure>\d+\.\d*)\s*Atm\.",
        description="The temperature and pressure of the Gaussian calculation. link 716",
    )
    ROTATIONAL_TEMPERATURE = MolOPPattern(
        content_pattern=r"^\s*Rotational temperatures \(Kelvin\)(?P<a>\s*\d+\.\d*)"
        r"(?P<b>\s*\d+\.\d*)(?P<c>\s*\d+\.\d*)",
        description="The rotational temperatures of the Gaussian calculation. link 716",
    )
    ROTATIONAL_CONST_IN_FREQUENCY_ANALYSIS = MolOPPattern(
        content_pattern=r"^\s*Rotational constants \(GHZ\):\s*"
        rf"(?P<a>{_ROTATIONAL_CONSTANT})\s+"
        rf"(?P<b>{_ROTATIONAL_CONSTANT})\s+"
        rf"(?P<c>{_ROTATIONAL_CONSTANT})",
        description="The rotational constants of the Gaussian calculation. link 716",
    )
    VIBRATIONAL_TEMPERATURE = MolOPPattern(
        content_pattern=r"^\s*(?:Vibrational temperatures:|\(Kelvin\)|)"
        r"(?P<values>(?:\s*\d+\.\d*){1,5})\n",
        content_repeat=0,
        description="The vibrational temperatures of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_PART = MolOPPattern(
        start_pattern="Zero-point correction",
        start_regex=False,
        end_pattern=r"^\s*Rotational(?:\s*-*\d+.\d+D[+-]\d+)\s*(?:-?\d+\.\d*)\s*(?:-?\d+\.\d*)",
        description="The thermochemistry part of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_CORRECTION = MolOPPattern(
        content_pattern=r"^\s*(?P<kind>Zero-point|Thermal) correction"
        r"(?P<suffix>.*)=\s*(?P<value>-?\d+\.\d*)",
        content_repeat=0,
        description="The thermochemistry corrections of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_SUM = MolOPPattern(
        content_pattern=r"^\s*Sum of electronic and (?P<term>thermal Free Energies|"
        r"thermal Enthalpies|thermal Energies|zero-point Energies)=\s*(?P<value>-?\d+.\d*)",
        content_repeat=0,
        description="The thermochemistry sums of the Gaussian calculation. link 716",
    )
    THERMOCHEMISTRY_CV_S = MolOPPattern(
        start_pattern="E (Thermal)             CV                S",
        start_regex=False,
        content_pattern=r"^\s*Total(?P<thermal_energy>\s*-?\d+\.\d*)"
        r"(?P<cv>\s*-?\d+\.\d*)(?P<entropy>\s*-?\d+\.\d*)",
        description="The thermochemistry CV and S of the Gaussian calculation. link 716",
    )
    MOLECULAR_MASS = MolOPPattern(
        content_pattern=r"^\s*Molecular mass:\s*(?P<mass>\d+\.\d+)\s*amu\.",
        description="The molecular mass in the thermochemistry section. link 716",
    )
    MOMENTS_OF_INERTIA = MolOPPattern(
        start_pattern="Principal axes and moments of inertia in atomic units:",
        start_regex=False,
        content_pattern=r"^\s*Eigenvalues --\s*(?P<values>.+)",
        description="The principal moments of inertia in atomic units. link 716",
    )
    ROTATIONAL_SYMMETRY_NUMBER = MolOPPattern(
        content_pattern=r"^\s*Rotational symmetry number\s+(?P<number>\d+)\.",
        description="The rotational symmetry number in the thermochemistry section. link 716",
    )
    FORCES_IN_CARTESIAN = MolOPPattern(
        start_pattern="Center     Atomic                   Forces (Hartrees/Bohr)",
        start_regex=False,
        end_pattern=r"^\s*-{3,}",
        end_offset=1,
        content_pattern=r"^\s*\d+\s+\d+(?P<x>\s*-?\d+\.\d*)"
        r"(?P<y>\s*-?\d+\.\d*)(?P<z>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The forces in Cartesian coordinates of the Gaussian calculation. link 716",
    )
    HESSIAN_IN_CARTESIAN = MolOPPattern(
        start_pattern="Force constants in Cartesian coordinates",
        start_regex=False,
        end_pattern=r"^\s*FormGI is forming",
        content_pattern=r"^\s*(?P<row>\d+)(?P<values>(?:\s*-?\d+.\d+D[+-]\d+){1,5})",
        content_repeat=0,
        description="The Hessian in Cartesian coordinates of the Gaussian calculation. link 716",
    )
    BERNY_STATE_MAJOR_PART = MolOPPattern(
        start_pattern="Item               Value     Threshold  Converged?",
        start_regex=False,
        end_pattern=r"(?r)^\s*^\s*Leave Link  103.*MaxMem=\s*(?:\d+)\s*cpu:\s*(?:\d+\.\d*)",
        description="The Berny optimization state. Patterns showed when #p used. link 103",
    )
    BERNY_STATE_BACKUP_PART = MolOPPattern(
        start_pattern="Item               Value     Threshold  Converged?",
        start_regex=False,
        end_pattern=r"^\s*(?:Grad){3,}",
        description="The Berny optimization state. link 103",
    )
    BERNY_STATE = MolOPPattern(
        content_pattern=r"(?r)^\s*(?P<label>Maximum\s+Force|RMS\s+Force|"
        r"Maximum\s+Displacement|RMS\s+Displacement)\s+(?P<value>\d+\.\d*)\s+"
        r"(?P<threshold>\d+\.\d*)\s+(?P<converged>NO|YES)",
        content_repeat=-4,
        description="The Berny optimization state. link 103",
    )
    ENERGY_CHANGE = MolOPPattern(
        content_pattern=r"(?r)^\s*^\s*Predicted change in Energy=(?P<value>\s*-?\d+.\d+D[+-]\d+)",
        description="The energy change of the Berny optimization state. link 103",
    )
    BERNY_CONCLUSION = MolOPPattern(
        content_pattern=r"(?r)^\s*^\s*Optimization completed\.",
        description="The Berny optimization conclusion. link 103",
    )
    ELECTRIC_DIPOLE_PART = MolOPPattern(
        start_pattern="Electric dipole moment (input orientation):",
        start_regex=False,
        end_pattern=r"^\s*-{3,}",
        description="The electric dipole moment part of the Gaussian calculation. link 9999",
    )
    ELECTRIC_DIPOLE_MOMENT = MolOPPattern(
        start_pattern="Electric dipole moment (input orientation):",
        start_regex=False,
        content_pattern=r"^\s*(?P<component>Tot|x|y|z)(?P<eigenvalue>\s*-?\d+\.\d*D[+-]\d+)"
        r"(?P<x>\s*-?\d+\.\d*D[+-]\d+)(?P<y>\s*-?\d+\.\d*D[+-]\d+)",
        content_repeat=4,
        description="The electric dipole moment of the Gaussian calculation. link 9999",
    )
    DIPOLE_POLARIZABILITY = MolOPPattern(
        start_pattern="Dipole polarizability, Alpha (input orientation)",
        start_regex=False,
        end_pattern=r"^\s*-{3,}",
        content_pattern=r"^\s*(?P<component>iso|aniso|xx|yx|yy|zx|zy|zz)"
        r"(?P<frequency_0>\s*-?\d+\.\d*D[+-]\d+)"
        r"(?P<frequency_1>\s*-?\d+\.\d*D[+-]\d+)"
        r"(?P<frequency_2>\s*-?\d+\.\d*D[+-]\d+)",
        content_repeat=8,
        description="The dipole polarizability of the Gaussian calculation. link 9999",
    )
    ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"^\s*1[\\|]1[\\|]GINC",
        end_pattern=r"(?:\\@\n|@\n)",
        description="The archive tail of the Gaussian calculation. link 9999",
    )
    JOB_TYPE_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        start_regex=False,
        start_offset=2,
        end_pattern="\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\(?P<value>.*)\\",
        description="The job type in the 4th block of the archive tail of the Gaussian calculation. link 9999",
    )
    FUNCTIONAL_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        start_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\(?P<value>.*)\\",
        description="The functional in the 5th block of the archive tail of the Gaussian calculation. link 9999",
    )
    BASIS_SET_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        start_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\(?P<value>.*)\\",
        description="The basis set in the 6th block of the archive tail of the Gaussian calculation. link 9999",
    )
    KEYWORDS_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        start_regex=False,
        end_pattern="\\\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\\\(?P<value>.*)\\\\",
        description="The keywords in the 12th block of the archive tail of the Gaussian calculation. link 9999",
    )
    TITLE_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        start_regex=False,
        end_pattern="\\\\",
        end_regex=False,
        end_offset=1,
        content_pattern=r"\\\\(?P<value>.*)\\\\",
        description="The title in the 14th block of the archive tail of the Gaussian calculation. link 9999",
    )
    CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        start_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=2,
        content_pattern=r"\\(?P<charge>-?\d+),(?P<multiplicity>-?\d+)\\",
        description="The charge and spin multiplicity in the 16th block of the archive tail of the Gaussian calculation. link 9999",
    )
    COORS_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\",
        start_regex=False,
        end_pattern="\\\\",
        end_regex=False,
        content_pattern=r"\\(?P<symbol>[A-Z][a-z]?),(?P<x>-?\d+\.\d*),"
        r"(?P<y>-?\d+\.\d*),(?P<z>-?\d+\.\d*)",
        content_repeat=0,
        description="The coordinates in the 17th block of the archive tail of the Gaussian calculation. link 9999",
    )
    VERSION_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern="\\\\",
        start_regex=False,
        end_pattern="\\",
        end_regex=False,
        end_offset=2,
        content_pattern=r"\\\\Version=(?P<value>.*)\\",
        description="The version in the archive tail of the Gaussian calculation. link 9999",
    )
    ENERGIES_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"(?P<method>HF|MP2|MP3|MP4[SDTQ]*|CCSD[\(\)T]*)="
        r"(?P<energy>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The energies in the archive tail of the Gaussian calculation. link 9999",
    )
    THERMOCHEMISTRY_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"(?P<term>ZeroPoint|Thermal|ETot|HTot|GTot)="
        r"(?P<value>\s*-?\d+\.\d*)",
        content_repeat=0,
        description="The thermochemistry in the archive tail of the Gaussian calculation. link 9999",
    )
    DIPOLE_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Dipole\s*=\s*(?P<x>-?\d+\.\d*)\s*,\s*"
        r"(?P<y>-?\d+\.\d*)\s*,\s*(?P<z>-?\d+\.\d*)",
        description="The dipole in the archive tail of the Gaussian calculation. link 9999",
    )
    POLAR_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Polar\s*=\s*(?P<xx>-?\d+\.\d*)\s*,\s*"
        r"(?P<xy>-?\d+\.\d*)\s*,\s*(?P<yy>-?\d+\.\d*)\s*,\s*"
        r"(?P<xz>-?\d+\.\d*)\s*,\s*(?P<yz>-?\d+\.\d*)\s*,\s*"
        r"(?P<zz>-?\d+\.\d*)",
        description="The polarizability in the archive tail of the Gaussian calculation. link 9999",
    )
    QUADRUPOLE_IN_ARCHIVE_TAIL = MolOPPattern(
        content_pattern=r"Quadrupole\s*=\s*(?P<xx>-?\d+\.\d*)\s*,\s*"
        r"(?P<yy>-?\d+\.\d*)\s*,\s*(?P<zz>-?\d+\.\d*)\s*,\s*"
        r"(?P<xy>-?\d+\.\d*)\s*,\s*(?P<xz>-?\d+\.\d*)\s*,\s*"
        r"(?P<yz>-?\d+\.\d*)",
        description="The quadrupole in the archive tail of the Gaussian calculation. link 9999",
    )
    HESSIAN_IN_ARCHIVE_TAIL = MolOPPattern(
        start_pattern=r"NImag=\d+\\\\",
        content_pattern=r"(?P<value>-?\d+\.\d*)",
        end_pattern=r"\\\\",
        end_offset=1,
        content_repeat=0,
        description="The Hessian in the archive tail of the Gaussian calculation. link 9999",
    )
    JOB_TIME = MolOPPattern(
        content_pattern=r"(?r)^\s*(?P<label>Job cpu time|Elapsed time):\s*"
        r"(?P<days>\d+)\s*days\s*(?P<hours>\d+)\s*hours\s*"
        r"(?P<minutes>\d+)\s*minutes\s*(?P<seconds>\d+\.\d*)\s*seconds",
        content_repeat=0,
        description="The job cpu or elapsed time of the Gaussian calculation. link 9999",
    )
    TERMINATION_STATUS = MolOPPattern(
        content_pattern=r"(?r)^\s*(?P<status>Normal|Error) termination (?P<message>.+)\s*"
        r"(?P<weekday>Mon|Tue|Wed|Thu|Fri|Sat|Sun)\s*"
        r"(?P<month>Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s*"
        r"(?P<day>\d+)\s*(?P<hour>\d+)\s*:\s*(?P<minute>\d+)\s*:\s*"
        r"(?P<second>\d+)\s*(?P<year>\d+)",
        description="The termination status of the Gaussian calculation. link 9999",
    )
    PROCEDURE_TIME = MolOPPattern(
        content_pattern=r"MaxMem=\s*(?P<max_mem>\d+)\s*cpu:\s*(?P<cpu>\d+.\d+)\s*"
        r"elap:\s*(?P<elapsed>\d+.\d+)",
        content_repeat=0,
        description="The procedure time of the end of link. link any",
    )


g16_log_patterns = G16LogPatterns()
