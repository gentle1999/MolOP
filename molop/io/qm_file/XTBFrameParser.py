"""
Author: TMJ
Date: 2024-06-21 11:03:55
LastEditors: TMJ
LastEditTime: 2024-06-25 20:53:21
Description: 请填写简介
"""

import os
import re

import numpy as np
from openbabel import pybel
from packaging.version import Version
from pydantic import Field, computed_field
from rdkit import Chem

from molop.io.bases.BaseMolFrameParser import BaseQMMolFrameParser
from molop.io.bases.DataClasses import (
    BondOrders,
    ChargeSpinPopulations,
    Energies,
    MolecularOrbitals,
    Polarizability,
    SinglePointProperties,
    ThermalEnergies,
)
from molop.logger.logger import moloplogger
from molop.unit import atom_ureg
from molop.utils.xtbpatterns import xtboutpatterns


class XTBFrameParser(BaseQMMolFrameParser):
    _block_type = "XTB OUT"
    qm_software: str = Field(default="xTB")
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    n_atom: int = Field(default=0, exclude=True, repr=False)

    def _parse(self):
        self.__block = self.frame_content
        self._parse_time()
        self._parse_coords()
        if self.only_extract_structure:
            return
        self._parse_state()
        self._parse_energy()
        self._parse_charges()
        self._parse_orbitals()
        self._parse_bond_order()
        self._parse_polarizability()
        self._parse_thermal()
        self._parse_rotation()
        self.single_point_properties = SinglePointProperties.model_validate(
            {
                **self._parse_vipea().model_dump(exclude_unset=True),
                **self._parse_gei().model_dump(exclude_unset=True),
                **self._parse_fukui_index().model_dump(exclude_unset=True),
            }
        )

    def _parse_time(self):
        if time_match := xtboutpatterns["time"].findall(self.__block):
            self.running_time = (
                sum(
                    float(day) * 24 * 60 * 60
                    + float(hour) * 60 * 60
                    + float(minute) * 60
                    + float(second)
                    for t, day, hour, minute, second in time_match
                )
                * atom_ureg.second
            )

    def _parse_coords_attached(self):
        attached_coords_path = re.search(
            xtboutpatterns["attached_coords"], self.__block
        )
        if attached_coords_path:
            coords_path = attached_coords_path.group(1)
            base_path = os.path.basename(coords_path)
            coords_type = coords_path.split(".")[-1]
            attach_path = None

            for file_path in [
                os.path.join(self.file_dir_path, base_path),
                os.path.join(self.file_dir_path, coords_path),
                os.path.join(self.file_dir_path, f"{self.pure_filename}.{coords_type}"),
                coords_path,
            ]:
                if os.path.exists(file_path):
                    attach_path = file_path
                    break

            if attach_path is None:
                moloplogger.error(
                    f"No attached coordinates found. Check your path if space in it. {coords_path}"
                )
                return False

            ofile = pybel.readfile(coords_type, attach_path)
            omol = next(ofile)
            self.coords = (
                np.array([atom.coords for atom in omol.atoms]) * atom_ureg.angstrom
            )
            self.atoms = [atom.atomicnum for atom in omol.atoms]
        else:
            moloplogger.error(
                "No attached coordinates found. Check your path if space in it."
            )
            return False
        return True

    def _parse_coords(self):
        coords_match = xtboutpatterns["coords_start"].search(self.__block)
        if not coords_match:
            if not self._parse_coords_attached():
                if "--hess" in self.keywords:
                    moloplogger.error(
                        "No coordinates found. If you are using `--hess`, use `--ohess` instead. See https://xtb-docs.readthedocs.io/en/latest/hessian.html for more details."
                    )
                raise RuntimeError("No coordinates found.")
        else:
            coords_type_match = xtboutpatterns["coords_type"].search(self.__block)
            if coords_type_match:
                coords_type = coords_type_match.group(1)
            coords_end = re.search(xtboutpatterns["coords_end"], self.__block)
            coords_block = self.__block[coords_match.start() : coords_end.end()]
            if coords_type == "xyz":
                if self.qm_software_version is None:
                    version = Version("0.0.0")
                else:
                    version = Version(self.qm_software_version.split()[1])
                if version <= Version("6.2.2"):
                    self._parse_coords_old(coords_block)
                else:
                    self._parse_coords_new(coords_block)
            else:
                omol = pybel.readstring(
                    coords_type, "\n".join(coords_block.split("\n")[2:])
                )
                self.coords = (
                    np.array([atom.coords for atom in omol.atoms]) * atom_ureg.angstrom
                )
                self.atoms = [atom.atomicnum for atom in omol.atoms]

        self.n_atom = len(self.atoms)

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
            self.atoms.append(Chem.Atom(atom).GetAtomicNum())
        self.coords = np.array(temp_coords) * atom_ureg.bohr

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
            self.atoms.append(Chem.Atom(atom).GetAtomicNum())
        self.coords = np.array(temp_coords) * atom_ureg.angstrom

    def _parse_state(self):
        if "opt" in self.keywords:
            if xtboutpatterns["geometric_optimization_state"].search(self.__block):
                self.geometry_optimization_status.geometry_optimized = True
            if match := xtboutpatterns["energy_convergence"].search(self.__block):
                self.geometry_optimization_status.energy_change_threshold = float(
                    match.group(1)
                )
            if match := xtboutpatterns["grad_convergence"].search(self.__block):
                self.geometry_optimization_status.max_force_threshold = float(
                    match.group(1)
                )
            if match := xtboutpatterns["energy_change"].findall(self.__block):
                self.geometry_optimization_status.energy_change = abs(float(match[-1]))
            if match := xtboutpatterns["grad_norm"].findall(self.__block):
                self.geometry_optimization_status.max_force = abs(float(match[-1]))
        if self.energies.total_energy is not None:
            self.status.scf_converged = True
        if xtboutpatterns["success_tag"].search(self.__block):
            self.status.normal_terminated = True

    @computed_field
    @property
    def is_error(self) -> bool:
        if not self.status.normal_terminated:
            return True
        if self.energies.total_energy is None:
            return True
        if (
            "opt" in self.keywords
            and self.geometry_optimization_status.geometry_optimized == False
        ):
            return True
        return False

    @property
    def is_optimized(self) -> bool:
        return self.status.normal_terminated

    def _parse_energy(self):
        energies = {}
        block = self.__block
        e_match = xtboutpatterns["total_energy"].search(block)
        while e_match:
            energies["scf_energy"] = (
                round(float(e_match.group(1)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            block = block[e_match.end() :]
            e_match = xtboutpatterns["total_energy"].search(block)
        block = self.__block
        e_match = xtboutpatterns["total_E"].search(block)
        while e_match:
            energies["scf_energy"] = (
                round(float(e_match.group(2)), 6)
                * atom_ureg.hartree
                / atom_ureg.particle
            )
            block = block[e_match.end() :]
            e_match = xtboutpatterns["total_E"].search(block)
        if len(energies):
            self.energies = Energies.model_validate(energies)

    def _parse_charges(self):
        charges = {}
        charges_match = xtboutpatterns["gfn1_charges_start"].search(self.__block)
        if charges_match:
            mulliken_charges = []
            cm5_charges = []
            for row in self.__block[charges_match.end() :].splitlines()[: self.n_atom]:
                mulliken_charges.append(float(row.split()[1]))
                cm5_charges.append(float(row.split()[2]))
            charges["mulliken_charges"] = mulliken_charges
            charges["hirshfeld_q_cm5"] = cm5_charges
        charges_match = xtboutpatterns["gfn2_charges_start"].search(self.__block)
        if charges_match:
            mulliken_charges = []
            for row in self.__block[charges_match.end() :].splitlines()[: self.n_atom]:
                mulliken_charges.append(float(row.split()[4]))
            charges["mulliken_charges"] = mulliken_charges
        if len(charges):
            self.charge_spin_populations = ChargeSpinPopulations.model_validate(charges)

    def _parse_orbitals(self):
        orbitals_start = xtboutpatterns["orbitals_start"].search(self.__block)
        if orbitals_start:
            orbitals_end = xtboutpatterns["orbitals_end"].search(
                self.__block[orbitals_start.end() :]
            )
        if orbitals_start and orbitals_end:
            block = self.__block[
                orbitals_start.end() : orbitals_start.end() + orbitals_end.start()
            ]
            orbital_energies = []
            alpha_orbital_occ = []
            beta_orbital_occ = []
            orbital_matches = xtboutpatterns["orbitals_match"].findall(block)
            has_homo = False
            for idx, occ, energy_Eh, energr_eV, tag in orbital_matches:
                if idx == "...":
                    continue
                if occ in (" ", "..."):
                    occ_value = 0
                else:
                    occ_value = int(float(occ))
                if int(idx) > len(orbital_energies) + 1:
                    for i in range(int(idx) - 1 - len(orbital_energies)):
                        orbital_energies.append(float("-inf"))
                        alpha_orbital_occ.append(not has_homo)
                        beta_orbital_occ.append(not has_homo)
                orbital_energies.append(float(energy_Eh))
                alpha_orbital_occ.append(occ_value > 0)
                beta_orbital_occ.append(occ_value > 1)
                if tag == "(HOMO)":
                    has_homo = True
            orbital_energies = (
                np.array(orbital_energies) * atom_ureg.hartree / atom_ureg.particle
            )
            self.molecular_orbitals = MolecularOrbitals(
                alpha_energies=orbital_energies,
                beta_energies=orbital_energies,
                alpha_occupancies=alpha_orbital_occ,
                beta_occupancies=beta_orbital_occ,
            )

    def _parse_bond_order(self):
        bond_order = {}
        start_matches = xtboutpatterns[f"wiberg_start"].search(self.__block)
        end_matches = xtboutpatterns[f"wiberg_end"].search(self.__block)
        result = np.zeros((self.n_atom, self.n_atom))
        if start_matches:
            temp_block = self.__block[start_matches.start() : end_matches.end()]
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
            bond_order["wiberg_bond_order"] = result
            self.bond_orders = BondOrders.model_validate(bond_order)

    def _parse_polarizability(self):
        polar = {}
        dipole_start = xtboutpatterns["dipole_start"].search(self.__block)
        dipole_end = xtboutpatterns["dipole_end"].search(self.__block)
        if dipole_start and dipole_end:
            temp_block = self.__block[dipole_start.start() : dipole_end.end()]
            values = list(map(float, xtboutpatterns["value"].findall(temp_block)))
            polar["dipole"] = (
                (np.array(values[:3]) + np.array(values[3:6])).astype(np.float32)
                * atom_ureg.atomic_unit_of_current
                * atom_ureg.atomic_unit_of_time
                * atom_ureg.bohr
            )
            self.polarizability = Polarizability.model_validate(polar)

    def _parse_thermal(self):
        thermal = {}
        zpc = xtboutpatterns["zp correct"].search(self.__block)
        if zpc:
            thermal["ZPVE"] = (
                float(zpc.group(1)) * atom_ureg.hartree / atom_ureg.particle
            )
        hc = xtboutpatterns["H correct"].search(self.__block)
        if hc:
            thermal["TCH"] = float(hc.group(1)) * atom_ureg.hartree / atom_ureg.particle
        gc = xtboutpatterns["G correct"].search(self.__block)
        if gc:
            thermal["TCG"] = float(gc.group(1)) * atom_ureg.hartree / atom_ureg.particle
        sh = xtboutpatterns["sum H"].search(self.__block)
        if sh:
            thermal["H_T"] = float(sh.group(1)) * atom_ureg.hartree / atom_ureg.particle
        sg = xtboutpatterns["sum G"].search(self.__block)
        if sg:
            thermal["G_T"] = float(sg.group(1)) * atom_ureg.hartree / atom_ureg.particle
        thermal_Cv_S_match = xtboutpatterns["thermal_Cv_S_start"].search(self.__block)
        if thermal_Cv_S_match:
            block = self.__block[thermal_Cv_S_match.end() :]
            for line in block.splitlines():
                if "TOT" in line:
                    thermal["C_V"] = (
                        float(line.split()[-3])
                        * atom_ureg.calorie
                        / atom_ureg.mol
                        / atom_ureg.kelvin
                    )
                    thermal["S"] = (
                        float(line.split()[-2])
                        * atom_ureg.calorie
                        / atom_ureg.mol
                        / atom_ureg.kelvin
                    )
                    break
        if len(thermal):
            self.thermal_energies = ThermalEnergies.model_validate(thermal)

    def _parse_rotation(self):
        rot = xtboutpatterns["rotation_consts"].search(self.__block)
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
            )

    def _parse_freqs(self):
        pass

    def _parse_vipea(self):
        single_property = {}
        vip_match = xtboutpatterns["VIP"].search(self.__block)
        if vip_match:
            single_property["vip"] = (
                float(vip_match.group(1)) * atom_ureg.eV / atom_ureg.particle
            )
        vea_match = xtboutpatterns["VEA"].search(self.__block)
        if vea_match:
            single_property["vea"] = (
                float(vea_match.group(1)) * atom_ureg.eV / atom_ureg.particle
            )
        return SinglePointProperties.model_validate(single_property)

    def _parse_gei(self):
        single_property = {}
        gei_match = xtboutpatterns["GEI"].search(self.__block)
        if gei_match:
            single_property["gei"] = (
                float(gei_match.group(1)) * atom_ureg.eV / atom_ureg.particle
            )
        return SinglePointProperties.model_validate(single_property)

    def _parse_fukui_index(self):
        single_property = {}
        fukui_match = xtboutpatterns["fukui_start"].search(self.__block)
        fukui_plus = []
        fukui_minus = []
        fukui_zero = []
        if fukui_match:
            for row in self.__block[fukui_match.end() :].splitlines()[: self.n_atom]:
                fukui_plus.append(float(row.split()[1]))
                fukui_minus.append(float(row.split()[2]))
                fukui_zero.append(float(row.split()[3]))
            single_property["fukui_positive"] = fukui_plus
            single_property["fukui_negative"] = fukui_minus
            single_property["fukui_zero"] = fukui_zero
        return SinglePointProperties.model_validate(single_property)
