import re
import os
from typing import Literal
from openbabel import pybel

import numpy as np
from packaging.version import Version

from molop.io.bases.molblock_base import QMBaseBlockParser
from molop.logger.logger import logger
from molop.unit import atom_ureg
from molop.utils.xtbpatterns import xtboutpatterns


class XTBOUTBlockParser(QMBaseBlockParser):
    """
    Parser for xTB out Blocks.
    """

    _block_type = "XTB OUT"

    def __init__(
        self,
        block: str,
        charge=0,
        multiplicity=1,
        basis: str = None,
        functional: str = None,
        solvent_model: str = None,
        solvent: str = None,
        temperature: float = 298.15,
        file_path="",
        version=None,
        parameter_comment=None,
        only_extract_structure=False,
    ):
        super().__init__(block, only_extract_structure)
        self.qm_software = "xTB"
        self._file_path = file_path
        self._charge = charge
        self._multiplicity = multiplicity
        self.version = version
        self.basis = basis
        self.functional = functional
        self.solvent_model = solvent_model
        self.solvent = solvent
        self.temperature = temperature
        self.parameter_comment = parameter_comment
        self._parse_coords()
        if not self._only_extract_structure:
            self._parse()

    def _parse(self):
        self._parse_state()
        self._parse_energy()
        self._parse_charges()
        self._parse_orbitals()
        self._parse_bond_order()
        self._parse_dipole()
        self._parse_thermo()
        self._parse_rotation()
        self._parse_vipea()
        self._parse_gei()
        self._parse_fukui_index()

    def _parse_coords_attached(self):
        attached_coords_path = re.search(xtboutpatterns["attached_coords"], self._block)
        if attached_coords_path:
            coords_path = attached_coords_path.group(1)
            base_path = os.path.basename(coords_path)
            coords_type = coords_path.split(".")[-1]
            attach_path = None

            for file_path in [
                os.path.join(self.file_dir_path, base_path),
                os.path.join(self.file_dir_path, coords_path),
                os.path.join(self.file_dir_path, f"{self.pure_filename}.{coords_type}"),
            ]:
                if os.path.exists(file_path):
                    attach_path = file_path
                    break

            if attach_path is None:
                logger.error(
                    f"No attached coordinates found. Check your path if space in it. {coords_path}"
                )
                return False

            ofile = pybel.readfile(coords_type, attach_path)
            omol = next(ofile)
            self._coords = (
                np.array([atom.coords for atom in omol.atoms]) * atom_ureg.angstrom
            )
            self._atoms = [atom.atomicnum for atom in omol.atoms]
        else:
            logger.error(
                "No attached coordinates found. Check your path if space in it."
            )
            return False

        return True

    def _parse_coords(self):
        coords_match = xtboutpatterns["coords_start"].search(self._block)
        if not coords_match:
            if not self._parse_coords_attached():
                if "--hess" in self.parameter_comment:
                    logger.error(
                        "No coordinates found. If you are using `--hess`, use `--ohess` instead. See https://xtb-docs.readthedocs.io/en/latest/hessian.html for more details."
                    )
                raise RuntimeError("No coordinates found.")
        else:
            coords_type_match = xtboutpatterns["coords_type"].search(self._block)
            if coords_type_match:
                coords_type = coords_type_match.group(1)
            coords_end = re.search(xtboutpatterns["coords_end"], self._block)
            coords_block = self._block[coords_match.start() : coords_end.end()]
            if coords_type == "xyz":
                if self.version is None:
                    version = Version("0.0.0")
                else:
                    version = Version(self.version.split()[1])
                if version <= Version("6.2.2"):
                    self._parse_coords_old(coords_block)
                else:
                    self._parse_coords_new(coords_block)
            else:
                omol = pybel.readstring(
                    coords_type, "\n".join(coords_block.split("\n")[2:])
                )
                self._coords = (
                    np.array([atom.coords for atom in omol.atoms]) * atom_ureg.angstrom
                )
                self._atoms = [atom.atomicnum for atom in omol.atoms]

        self.__n_atom = len(self.atoms)

    def _parse_coords_old(self, block: str):
        """
        e.g.

        ================
         final structure:
        ================
        $coord
            2.52072290250473   -0.04782551206377   -0.50388676977877      C
                    .                 .                    .              .
        """
        coords = xtboutpatterns["coords_old"].findall(block)
        temp_coords = []
        for x, y, z, atom in coords:
            temp_coords.append((float(x), float(y), float(z)))
            self._atoms.append(atom)
        self._coords = (np.array(temp_coords) * atom_ureg.bohr).to("angstrom")

    def _parse_coords_new(self, block: str):
        """
        e.g.

        ================
         final structure:
        ================
        5
         xtb: 6.2.3 (830e466)
        Cl        1.62694523673790    0.09780349799138   -0.02455489507427
        C        -0.15839164427314   -0.00942638308615    0.00237760557913
        H        -0.46867957388620   -0.59222865914178   -0.85786049981721
        H        -0.44751262498645   -0.49575975568264    0.92748366742968
        H        -0.55236139359212    0.99971129991918   -0.04744587811734
        """
        coords = xtboutpatterns["coords_new"].findall(block)
        temp_coords = []
        for atom, x, y, z in coords:
            temp_coords.append((float(x), float(y), float(z)))
            self._atoms.append(atom)
        self._coords = np.array(temp_coords) * atom_ureg.angstrom

    def _parse_state(self):
        if "opt" in self.parameter_comment:
            if xtboutpatterns["geometric_optimization_state"].search(self._block):
                self.status["geometric_optimization"] = True
            else:
                self.status["geometric_optimization"] = False
        if xtboutpatterns["success_tag"].search(self._block):
            self.status["finished run"] = True
        else:
            self.status["finished run"] = False

    def is_error(self) -> bool:
        if not self.status["finished run"]:
            return True
        if self.total_energy is None:
            return True
        if "geometric_optimization" in self.status:
            return self.status["geometric_optimization"] == False
        else:
            return True

    def is_optimized(self) -> bool:
        if "geometric_optimization" in self.status:
            return self.status["geometric_optimization"]
        return False

    def _parse_energy(self):
        block = self._block
        e_match = xtboutpatterns["total_energy"].search(block)
        while e_match:
            self._total_energy = (
                round(float(e_match.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            block = block[e_match.end() :]
            e_match = xtboutpatterns["total_energy"].search(block)
        block = self._block
        e_match = xtboutpatterns["total_E"].search(block)
        while e_match:
            self.scf_energy = (
                round(float(e_match.group(2)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            block = block[e_match.end() :]
            e_match = xtboutpatterns["total_E"].search(block)

    def _parse_charges(self):
        self.cm5_charges = []
        charges_match = xtboutpatterns["gfn1_charges_start"].search(self._block)
        if charges_match:
            mulliken_charges = []
            cm5_charges = []
            for row in self._block[charges_match.end() :].splitlines()[: self.__n_atom]:
                mulliken_charges.append(float(row.split()[1]))
                cm5_charges.append(float(row.split()[2]))
            self.mulliken_charges = mulliken_charges
            self.cm5_charges = cm5_charges
        charges_match = xtboutpatterns["gfn2_charges_start"].search(self._block)
        if charges_match:
            mulliken_charges = []
            for row in self._block[charges_match.end() :].splitlines()[: self.__n_atom]:
                mulliken_charges.append(float(row.split()[4]))
            self.mulliken_charges = mulliken_charges

    def _parse_orbitals(self):
        orbitals_start = xtboutpatterns["orbitals_start"].search(self._block)
        orbitals_end = xtboutpatterns["orbitals_end"].search(
            self._block[orbitals_start.end() :]
        )
        if orbitals_start and orbitals_end:
            block = self._block[
                orbitals_start.end() : orbitals_start.end() + orbitals_end.start()
            ]
            orbitals = []
            orbital_match = xtboutpatterns["orbitals_match"].search(block)
            while orbital_match:
                idx, occ, energy_Eh, energr_eV, tag = orbital_match.groups()
                if idx == "...":
                    block = block[orbital_match.end() :]
                    orbital_match = xtboutpatterns["orbitals_match"].search(block)
                    continue
                if int(idx) > len(orbitals) + 1:
                    for i in range(int(idx) - 1 - len(orbitals)):
                        orbitals.append(float("-inf"))
                orbitals.append(float(energy_Eh))
                if tag == "(HOMO)":
                    self.alpha_energy["homo"] = (
                        round(float(energy_Eh), 6)
                        * atom_ureg.hartree
                        / atom_ureg.particle
                    )
                if tag == "(LUMO)":
                    self.alpha_energy["lumo"] = (
                        round(float(energy_Eh), 6)
                        * atom_ureg.hartree
                        / atom_ureg.particle
                    )
                block = block[orbital_match.end() :]
                orbital_match = xtboutpatterns["orbitals_match"].search(block)
            self.alpha_energy["gap"] = (
                round((self.alpha_energy["lumo"] - self.alpha_energy["homo"]).m, 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            self.alpha_FMO_orbits = (
                np.array(orbitals) * atom_ureg.hartree / atom_ureg.particle
            )

    def _parse_bond_order(self):
        start_matches = xtboutpatterns[f"wiberg_start"].search(self._block)
        end_matches = xtboutpatterns[f"wiberg_end"].search(self._block)
        result = np.zeros((self.__n_atom, self.__n_atom))
        if start_matches:
            temp_block = self._block[start_matches.start() : end_matches.end()]
            for idx, sub_temp_block in enumerate(
                xtboutpatterns["index"].split(temp_block)[1:]
            ):
                bond = list(
                    map(
                        lambda x: int(x) - 1,
                        xtboutpatterns["bond"].findall(sub_temp_block),
                    )
                )
                value = list(
                    map(float, xtboutpatterns["value"].findall(sub_temp_block))
                )[1:]
                for i, v in zip(bond, value):
                    result[idx][i] = v
                    result[i][idx] = v
            self.wiberg_bond_order = result

    def _parse_dipole(self):
        dipole_start = xtboutpatterns["dipole_start"].search(self._block)
        dipole_end = xtboutpatterns["dipole_end"].search(self._block)
        if dipole_start and dipole_end:
            temp_block = self._block[dipole_start.start() : dipole_end.end()]
            values = list(map(float, xtboutpatterns["value"].findall(temp_block)))
            self._dipole = (np.array(values[:3]) + np.array(values[3:6])).astype(
                np.float32
            ) * atom_ureg.debye

    def _parse_thermo(self):
        zpc = xtboutpatterns["zp correct"].search(self._block)
        if zpc:
            self.thermal_energy["zpve"] = (
                float(zpc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        hc = xtboutpatterns["H correct"].search(self._block)
        if hc:
            self.thermal_energy["TCH"] = (
                float(hc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        gc = xtboutpatterns["G correct"].search(self._block)
        if gc:
            self.thermal_energy["TCG"] = (
                float(gc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        se = xtboutpatterns["sum E"].search(self._block)
        if se:
            self.thermal_energy["U_0"] = (
                float(se.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        sh = xtboutpatterns["sum H"].search(self._block)
        if sh:
            self.thermal_energy["H_T"] = (
                float(sh.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        sg = xtboutpatterns["sum G"].search(self._block)
        if sg:
            self.thermal_energy["G_T"] = (
                float(sg.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        thermal_Cv_S_match = xtboutpatterns["thermal_Cv_S_start"].search(self._block)
        if thermal_Cv_S_match:
            block = self._block[thermal_Cv_S_match.end() :]
            for line in block.splitlines():
                if "TOT" in line:
                    self.capacity = (
                        float(line.split()[-3])
                        * atom_ureg.calorie
                        / atom_ureg.mol
                        / atom_ureg.kelvin
                    )
                    self.entropy = (
                        float(line.split()[-2])
                        * atom_ureg.calorie
                        / atom_ureg.mol
                        / atom_ureg.kelvin
                    )
                    break

    def _parse_rotation(self):
        rot = xtboutpatterns["rotation_consts"].search(self._block)
        if rot:
            rots = []
            for i in range(1, 4):
                try:
                    temp_rot = float(rot.group(i))
                except:
                    temp_rot = 0
                rots.append(temp_rot)
            self.rotation_constants = (
                np.array(rots) * atom_ureg.cm_1 * 299792458 * atom_ureg.m / atom_ureg.s
            ).to("Ghertz")

    def _parse_freqs(self):
        pass

    def _parse_vipea(self):
        self.vip, self.vea = None, None
        vip_match = xtboutpatterns["VIP"].search(self._block)
        if vip_match:
            self.vip = float(vip_match.group(1)) * atom_ureg.eV / atom_ureg.particle
        vea_match = xtboutpatterns["VEA"].search(self._block)
        if vea_match:
            self.vea = float(vea_match.group(1)) * atom_ureg.eV / atom_ureg.particle

    def _parse_gei(self):
        self.gei = None
        gei_match = xtboutpatterns["GEI"].search(self._block)
        if gei_match:
            self.gei = float(gei_match.group(1)) * atom_ureg.eV / atom_ureg.particle

    def _parse_fukui_index(self):
        fukui_match = xtboutpatterns["fukui_start"].search(self._block)
        self.fukui_plus = []
        self.fukui_minus = []
        self.fukui_zero = []
        if fukui_match:
            for row in self._block[fukui_match.end() :].splitlines()[: self.__n_atom]:
                self.fukui_plus.append(float(row.split()[1]))
                self.fukui_minus.append(float(row.split()[2]))
                self.fukui_zero.append(float(row.split()[3]))
