"""
Author: TMJ
Date: 2023-12-16 21:29:31
LastEditors: TMJ
LastEditTime: 2024-02-07 15:41:40
Description: 请填写简介
"""
from openbabel import pybel
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
# ob_log_handler = pybel.ob.OBMessageHandler()
# ob_log_handler.SetOutputLevel(0)
pybel.ob.obErrorLog.StopLogging()

from molop.io import AutoParser
