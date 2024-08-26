'''
Author: TMJ
Date: 2024-01-09 10:38:19
LastEditors: TMJ
LastEditTime: 2024-08-26 19:48:15
Description: 请填写简介
'''
import logging

moloplogger = logging.getLogger("molop")
moloplogger.setLevel(logging.INFO)
file_handler = logging.FileHandler("molop.log")
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
sh_formatter = logging.Formatter("%(levelname)s - %(message)s")
stream_handler.setFormatter(sh_formatter)
moloplogger.addHandler(stream_handler)

moloplogger.addHandler(file_handler)