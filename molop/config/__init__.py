'''
Author: TMJ
Date: 2024-02-14 14:40:02
LastEditors: TMJ
LastEditTime: 2024-04-26 10:52:12
Description: 请填写简介
'''
import logging

from openbabel import pybel
from rdkit import RDLogger
from rdkit.Chem import AllChem

from molop.logger.logger import logger


class MolOPConfig:
    def __init__(self) -> None:
        self.show_progress_bar = True
        self.show_log = True
        self._sh = logging.StreamHandler()
        logger.addHandler(self._sh)
        RDLogger.DisableLog("rdApp.*")
        # ob_log_handler = pybel.ob.OBMessageHandler()
        # ob_log_handler.SetOutputLevel(0)
        pybel.ob.obErrorLog.StopLogging()
        self.max_jobs = 16
        self.morgan_fpgen = AllChem.GetMorganGenerator(radius=3, fpSize=1024, includeChirality=True)
        self.atompair_fpgen = AllChem.GetAtomPairGenerator(fpSize=1024, includeChirality=True)
        self.rdkit_fpgen = AllChem.GetRDKitFPGenerator(fpSize=1024)
        self.topological_torsion_fpgen = AllChem.GetTopologicalTorsionGenerator(fpSize=1024)
        self.max_structure_recovery_time = 10
        self.allow_spin_change = False

    def quiet(self):
        self.show_progress_bar = False
        self.show_log = False
        logger.removeHandler(self._sh)

    def verbose(self):
        self.show_progress_bar = True
        self.show_log = True
        logger.addHandler(self._sh)


molopconfig = MolOPConfig()
