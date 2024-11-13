"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2024-11-11 22:00:27
Description: 请填写简介
"""

import logging

moloplogger = logging.getLogger("molop")

"""
file_handler = logging.FileHandler("molop.log")
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
"""

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
sh_formatter = logging.Formatter("%(levelname)s - %(message)s")
stream_handler.setFormatter(sh_formatter)
moloplogger.addHandler(stream_handler)

moloplogger.setLevel(logging.INFO)
