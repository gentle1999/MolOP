'''
Author: TMJ
Date: 2024-01-09 10:38:19
LastEditors: TMJ
LastEditTime: 2024-08-30 11:32:56
Description: 请填写简介
'''
import logging

moloplogger = logging.getLogger("molop")

file_handler = logging.FileHandler("molop.log")
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
moloplogger.addHandler(file_handler)

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
sh_formatter = logging.Formatter("%(levelname)s - %(message)s")
stream_handler.setFormatter(sh_formatter)
moloplogger.addHandler(stream_handler)

"""debug_handler = logging.StreamHandler()
debug_handler.setLevel(logging.DEBUG)
debug_formatter = logging.Formatter("%(levelname)s - %(message)s")
debug_handler.setFormatter(debug_formatter)
moloplogger.addHandler(debug_handler)"""

moloplogger.setLevel(logging.INFO)