'''
Author: TMJ
Date: 2024-01-31 21:57:38
LastEditors: TMJ
LastEditTime: 2024-01-31 22:12:10
Description: 请填写简介
'''
from molop.logger.logger import logger
import logging

class MolOPConfig:
    def __init__(self) -> None:
        self.show_progress_bar = True
        self.show_log = True
        self._sh = logging.StreamHandler()
        logger.addHandler(self._sh)

    def quiet(self):
        self.show_progress_bar = False
        self.show_log = False
        logger.removeHandler(self._sh)

    def verbose(self):
        self.show_progress_bar = True
        self.show_log = True
        logger.addHandler(self._sh)


molopconfig = MolOPConfig()