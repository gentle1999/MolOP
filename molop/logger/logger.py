'''
Author: TMJ
Date: 2024-01-09 10:38:19
LastEditors: TMJ
LastEditTime: 2024-01-28 13:52:19
Description: 请填写简介
'''
import logging

logging.basicConfig(
    filename="molop.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

logger = logging.getLogger("molop")
sh = logging.StreamHandler()
logger.addHandler(sh)