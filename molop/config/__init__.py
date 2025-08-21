"""
Author: TMJ
Date: 2024-10-19 09:57:26
LastEditors: TMJ
LastEditTime: 2024-12-02 21:28:20
Description: 请填写简介
"""

import logging

from openbabel import pybel
from pydantic import BaseModel, ConfigDict, Field, ValidationError
from rdkit import RDLogger
from rdkit.Chem.rdFingerprintGenerator import (
    FingerprintGenerator64,
    GetAtomPairGenerator,
    GetMorganGenerator,
    GetRDKitFPGenerator,
    GetTopologicalTorsionGenerator,
)

moloplogger = logging.getLogger("moloplogger")
stream_handler = logging.StreamHandler()
moloplogger.addHandler(stream_handler)


class MolOPConfig(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    show_progress_bar: bool = Field(default=True, description="是否显示进度条")
    max_jobs: int = Field(default=16, description="最大作业数")
    morgan_fpgen: FingerprintGenerator64 = Field(
        default=GetMorganGenerator(radius=3, fpSize=1024, includeChirality=True),
        description="Morgan 指纹生成器",
    )
    atompair_fpgen: FingerprintGenerator64 = Field(
        default=GetAtomPairGenerator(fpSize=1024, includeChirality=True),
        description="AtomPair 指纹生成器",
    )
    rdkit_fpgen: FingerprintGenerator64 = Field(
        default=GetRDKitFPGenerator(fpSize=1024),
        description="RDKit 指纹生成器",
    )
    topological_torsion_fpgen: FingerprintGenerator64 = Field(
        default=GetTopologicalTorsionGenerator(fpSize=1024),
        description="TopologicalTorsion 指纹生成器",
    )
    strict_structure_recovery: bool = Field(
        default=False, description="是否进行严格的结构恢复"
    )
    max_structure_recovery_time: int = Field(default=10, description="最大结构恢复时间")
    allow_spin_change: bool = Field(default=False, description="是否允许自旋变化")
    force_unit_transform: bool = Field(default=False, description="是否强制单位转换")
    parallel_max_size: int = Field(default=8 * 1024**2, description="并行最大大小")

    def __init__(self, **data):
        super().__init__(**data)
        RDLogger.DisableLog("rdApp.*")  # type: ignore
        pybel.ob.obErrorLog.StopLogging()
        moloplogger.addHandler(stream_handler)

    def quiet(self):
        """
        Turn off progress bar and log messages to stdout.
        """
        if self.show_progress_bar:
            self.show_progress_bar = False
            moloplogger.removeHandler(stream_handler)

    def verbose(self):
        """
        Turn on progress bar and log messages to stdout.
        """
        if not self.show_progress_bar:
            self.show_progress_bar = True
            moloplogger.addHandler(stream_handler)

    def set_log_level(self, level: str):
        try:
            moloplogger.setLevel(level)
        except ValueError as e:
            logging.error(f"设置日志级别时出错: {e}")
            raise ValueError(f"无效的日志级别: {level}") from e


# 创建 MolOPConfig 实例
try:
    molopconfig = MolOPConfig()
except ValidationError as e:
    logging.error(f"配置验证失败: {e}")
    raise e
