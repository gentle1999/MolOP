"""
Author: TMJ
Date: 2024-02-14 14:40:02
LastEditors: TMJ
LastEditTime: 2024-03-20 15:04:21
Description: 请填写简介
"""
from molop.logger.logger import logger
import logging
from openbabel import pybel
from rdkit import RDLogger
from rdkit.Chem import AllChem


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
        self.fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    def quiet(self):
        self.show_progress_bar = False
        self.show_log = False
        logger.removeHandler(self._sh)

    def verbose(self):
        self.show_progress_bar = True
        self.show_log = True
        logger.addHandler(self._sh)


molopconfig = MolOPConfig()
