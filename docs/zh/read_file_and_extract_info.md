# 读取文件并提取信息

这里有一个如何读取文件并从中提取信息的例子。

假设在这样一个文件夹中有一些 g16log 文件：


```python
! ls -l ../../tests/test_files/g16log/
```

    total 33748
    -rw-r--r-- 1 cathayana cathayana  727184 Feb 20 11:47 11_Opt.log
    -rw-r--r-- 1 cathayana cathayana   53139 Feb 20 11:47 3_Sp.log
    -rw-r--r-- 1 cathayana cathayana 8008984 Feb 20 11:47 RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log
    -rw-r--r-- 1 cathayana cathayana 1424205 Feb 20 11:47 S_Ph_Ni_TS.log
    -rw-r--r-- 1 cathayana cathayana 1910081 Feb 20 11:47 TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log
    -rw-r--r-- 1 cathayana cathayana   65746 Feb 20 11:47 TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log
    -rw-r--r-- 1 cathayana cathayana 3896697 Feb 20 11:47 TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log
    -rw-r--r-- 1 cathayana cathayana 1307536 Feb 20 11:47 TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log
    -rw-r--r-- 1 cathayana cathayana  165117 Feb 20 11:47 dsgdb9nsd_000001-3-.log
    -rw-r--r-- 1 cathayana cathayana  823921 Feb 20 11:47 dsgdb9nsd_000107-3-.log
    -rw-r--r-- 1 cathayana cathayana  703600 Feb 20 11:47 dsgdb9nsd_000180-9-.log
    -rw-r--r-- 1 cathayana cathayana  385852 Feb 20 11:47 dsgdb9nsd_000484-1+.log
    -rw-r--r-- 1 cathayana cathayana  509600 Feb 20 11:47 dsgdb9nsd_000672-3+.log
    -rw-r--r-- 1 cathayana cathayana  583565 Feb 20 11:47 dsgdb9nsd_000696-4.log
    -rw-r--r-- 1 cathayana cathayana  380253 Feb 20 11:47 dsgdb9nsd_000763-2-.log
    -rw-r--r-- 1 cathayana cathayana  446864 Feb 20 11:47 dsgdb9nsd_000923-3+.log
    -rw-r--r-- 1 cathayana cathayana  523877 Feb 20 11:47 dsgdb9nsd_000955-3.log
    -rw-r--r-- 1 cathayana cathayana  383481 Feb 20 11:47 dsgdb9nsd_000958-3-.log
    -rw-r--r-- 1 cathayana cathayana  639463 Feb 20 11:47 dsgdb9nsd_001232-4-.log
    -rw-r--r-- 1 cathayana cathayana  535511 Feb 20 11:47 dsgdb9nsd_002924-8-.log
    -rw-r--r-- 1 cathayana cathayana  899961 Feb 20 11:47 dsgdb9nsd_003051-3.log
    -rw-r--r-- 1 cathayana cathayana  256920 Feb 20 11:47 dsgdb9nsd_003895.log
    -rw-r--r-- 1 cathayana cathayana  868095 Feb 20 11:47 dsgdb9nsd_004015-3-.log
    -rw-r--r-- 1 cathayana cathayana  442477 Feb 20 11:47 dsgdb9nsd_004478-6-.log
    -rw-r--r-- 1 cathayana cathayana  958696 Feb 20 11:47 dsgdb9nsd_004517-4-.log
    -rw-r--r-- 1 cathayana cathayana  672242 Feb 20 11:47 dsgdb9nsd_004669-4.log
    -rw-r--r-- 1 cathayana cathayana  453094 Feb 20 11:47 dsgdb9nsd_004738-2.log
    -rw-r--r-- 1 cathayana cathayana 1334751 Feb 20 11:47 dsgdb9nsd_006075rearrange.log
    -rw-r--r-- 1 cathayana cathayana  505226 Feb 20 11:47 dsgdb9nsd_009986.log
    -rw-r--r-- 1 cathayana cathayana  452946 Feb 20 11:47 dsgdb9nsd_130366.log
    -rw-r--r-- 1 cathayana cathayana  641086 Feb 20 11:47 dsgdb9nsd_131200-4-.log
    -rw-r--r-- 1 cathayana cathayana  499704 Feb 20 11:47 dsgdb9nsd_131200.log
    -rw-r--r-- 1 cathayana cathayana  911981 Feb 20 11:47 dsgdb9nsd_131941-4+.log
    -rw-r--r-- 1 cathayana cathayana  540614 Feb 20 11:47 dsgdb9nsd_132072.log
    -rw-r--r-- 1 cathayana cathayana  553016 Feb 20 11:47 dsgdb9nsd_133826.log
    -rw-r--r-- 1 cathayana cathayana  650970 Feb 20 11:47 dsgdb9nsd_133858.log
    -rw-r--r-- 1 cathayana cathayana  301634 Feb 20 11:47 molecule_0.log
    -rw-r--r-- 1 cathayana cathayana   22094 Feb 20 11:47 r1_C2H3N3O_sp_g16.log
    -rw-r--r-- 1 cathayana cathayana   38369 Feb 20 11:47 r1_C8H9N3O2_sp_g16.log


看起来很混乱？是的，它们是来自不同项目的混合物。MolOP 可以帮助你读取文件，并通过通用方法从中提取信息。


```python
from molop import AutoParser
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole  # for better drawing

IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 400, 400

files = AutoParser("../../tests/test_files/g16log/*.log")
```

    MolOP parsing with 28 jobs: 100%|██████████| 39/39 [00:03<00:00, 10.59it/s]
    0 files failed to parse, 39 successfully parsed


我们首先可以获得文件的摘要。摘要包含文件的结构信息（SMILES）和一些关键的 QM 信息。


```python
print(files.to_summary_df().to_markdown())
```

    |    | parser       | file_name                                                                    | file_path                                                                                                                 | file_format   |   charge |   multiplicity | SMILES                                                                                                 | status                                                                                                                                                                                                   |        ZPE |        TCE |        TCH |        TCG |   ZPE-Gas |     E-Gas |     H-Gas |     G-Gas |         sp |     HOMO |     LUMO |     GAP |   first freq | first freq tag   |   second freq | second freq tag   |     S**2 |        S |
    |---:|:-------------|:-----------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|:--------------|---------:|---------------:|:-------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------:|-----------:|-----------:|-----------:|----------:|----------:|----------:|----------:|-----------:|---------:|---------:|--------:|-------------:|:-----------------|--------------:|:------------------|---------:|---------:|
    |  0 | G16LOGParser | 11_Opt.log                                                                   | /home/tmj/proj/MolOP/tests/test_files/g16log/11_Opt.log                                                                   | .log          |        0 |              1 | C=C[C@H]1[C@@H]2[C@H](C[C@@]1(C)O)OC[C@@H]2C                                                           | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': True, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                    |   0.275585 |   0.287903 |   0.288847 |   0.238406 |   -579.65 |  -579.637 |  -579.636 |  -579.687 |  -579.925  | -0.2188  | -0.01452 | 0.20428 |       76.85  | False            |       94.3685 | False             | nan      | nan      |
    |  1 | G16LOGParser | 3_Sp.log                                                                     | /home/tmj/proj/MolOP/tests/test_files/g16log/3_Sp.log                                                                     | .log          |        0 |              1 | CCC[C@H]1CO[C@H]2C[C@@](C)(O)C[C@@H]12                                                                 | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -581.781  | -0.21593 |  0.02581 | 0.24174 |      nan     |                  |      nan      |                   | nan      | nan      |
    |  2 | G16LOGParser | RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log | /home/tmj/proj/MolOP/tests/test_files/g16log/RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log | .log          |       -1 |              1 | CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)CO2)[N-]1.C[CH-]C(=O)[N-]c1ccccc1.[Cu@TB7+5].[O-]c1ccccc1 | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': False, 'RMS     Displacement': True, 'termination': 'Normal', 'failure reason': '! Non-Optimized Parameters !', 'SCF Done': True} |   0.643965 |   0.682345 |   0.68329  |   0.572926 |  -1827.58 | -1827.54  | -1827.54  | -1827.65  | -1828.22   | -0.06449 |  0.0099  | 0.07439 |     -126.551 | True             |       14.3776 | False             | nan      | nan      |
    |  3 | G16LOGParser | S_Ph_Ni_TS.log                                                               | /home/tmj/proj/MolOP/tests/test_files/g16log/S_Ph_Ni_TS.log                                                               | .log          |        0 |              1 | [Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1.c1c[nH]cn1.c1c[nH]cn1                                    | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': True, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                    |   0.361561 |   0.38868  |   0.389624 |   0.301123 |  -1548.69 | -1548.66  | -1548.66  | -1548.75  | -1549.05   | -0.18821 | -0.11033 | 0.07788 |     -185.866 | True             |       15.8954 | False             |  -0      |  -0      |
    |  4 | G16LOGParser | TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log                                      | .log          |        0 |              1 | CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F                                                   | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': False, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                   |   0.251882 |   0.277379 |   0.278323 |   0.192899 |  -1287.11 | -1287.08  | -1287.08  | -1287.17  | -1287.36   | -0.25863 | -0.10047 | 0.15816 |     -300.831 | True             |       16.957  | False             | nan      | nan      |
    |  5 | G16LOGParser | TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log                                         | .log          |        0 |              1 | CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F                                                   | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     | -1288.85   | -0.26367 | -0.10703 | 0.15664 |      nan     |                  |      nan      |                   | nan      | nan      |
    |  6 | G16LOGParser | TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log                                    | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log                                    | .log          |        0 |              1 | C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1                                                      | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': False, 'RMS     Displacement': True, 'termination': 'Normal', 'failure reason': '! Non-Optimized Parameters !', 'SCF Done': True} |   0.349175 |   0.373165 |   0.374109 |   0.294821 |  -6028.57 | -6028.55  | -6028.55  | -6028.63  | -6028.92   | -0.20404 | -0.04169 | 0.16235 |      -56.863 | True             |       29.9487 | False             | nan      | nan      |
    |  7 | G16LOGParser | TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log                     | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log                     | .log          |        0 |              1 | CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(F)(F)F                                                   | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': True, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                    |   0.252491 |   0.277515 |   0.27846  |   0.195753 |  -1287.12 | -1287.09  | -1287.09  | -1287.17  | -1287.37   | -0.25447 | -0.11005 | 0.14442 |     -252.41  | True             |       25.8663 | False             | nan      | nan      |
    |  8 | G16LOGParser | dsgdb9nsd_000001-3-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000001-3-.log                                                      | .log          |       -1 |              1 | [CH3-]                                                                                                 | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |   -39.8575 |  0.06578 |  0.17809 | 0.11231 |      nan     |                  |      nan      |                   | nan      | nan      |
    |  9 | G16LOGParser | dsgdb9nsd_000107-3-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000107-3-.log                                                      | .log          |       -1 |              1 | [C-]#CCC#C                                                                                             | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -192.275  |  0.00382 |  0.13824 | 0.13442 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 10 | G16LOGParser | dsgdb9nsd_000180-9-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000180-9-.log                                                      | .log          |       -1 |              1 | CC/C(C)=N/[O-]                                                                                         | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -287.285  |  0.02146 |  0.11362 | 0.09216 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 11 | G16LOGParser | dsgdb9nsd_000484-1+.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000484-1+.log                                                      | .log          |        1 |              1 | [C+]#CC#CC#C                                                                                           | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -228.647  | -0.47931 | -0.43653 | 0.04278 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 12 | G16LOGParser | dsgdb9nsd_000672-3+.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000672-3+.log                                                      | .log          |        1 |              1 | [C+]1=C[C@H]2O[C@H]2C1                                                                                 | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -268.398  | -0.51954 | -0.32947 | 0.19007 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 13 | G16LOGParser | dsgdb9nsd_000696-4.log                                                       | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000696-4.log                                                       | .log          |        0 |              2 | C1[C@@H]2[C]3CN2[C@@H]13                                                                               | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -248.756  | -0.32784 | -0.18259 | 0.14525 |      nan     |                  |      nan      |                   |   0.7535 |   0.5018 |
    | 14 | G16LOGParser | dsgdb9nsd_000763-2-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000763-2-.log                                                      | .log          |       -1 |              1 | [H]/N=c1/[cH-]onn1                                                                                     | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -316.937  |  0.0093  |  0.13147 | 0.12217 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 15 | G16LOGParser | dsgdb9nsd_000923-3+.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000923-3+.log                                                      | .log          |        1 |              1 | O=C1C=C=NC=[NH+]1                                                                                      | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -338.692  | -0.4928  | -0.35012 | 0.14268 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 16 | G16LOGParser | dsgdb9nsd_000955-3.log                                                       | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000955-3.log                                                       | .log          |        0 |              2 | Oc1nc[c]cn1                                                                                            | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -338.98   | -0.28326 | -0.27967 | 0.00359 |      nan     |                  |      nan      |                   |   0.757  |   0.5035 |
    | 17 | G16LOGParser | dsgdb9nsd_000958-3-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_000958-3-.log                                                      | .log          |       -1 |              1 | Oc1n[c-]ncn1                                                                                           | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -355.094  | -0.00299 |  0.12011 | 0.1231  |      nan     |                  |      nan      |                   | nan      | nan      |
    | 18 | G16LOGParser | dsgdb9nsd_001232-4-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_001232-4-.log                                                      | .log          |       -1 |              1 | N/C=N/C=C(/[O-])O                                                                                      | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -377.445  |  0.00144 |  0.10845 | 0.10701 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 19 | G16LOGParser | dsgdb9nsd_002924-8-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_002924-8-.log                                                      | .log          |       -1 |              1 | COC1=NCC[N-]1                                                                                          | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -341.475  |  0.01    |  0.12531 | 0.11531 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 20 | G16LOGParser | dsgdb9nsd_003051-3.log                                                       | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_003051-3.log                                                       | .log          |        0 |              2 | O=CN[CH][C@H]1CO1                                                                                      | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -361.273  | -0.28613 | -0.19199 | 0.09414 |      nan     |                  |      nan      |                   |   0.7565 |   0.5032 |
    | 21 | G16LOGParser | dsgdb9nsd_003895.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_003895.log                                                         | .log          |        0 |              1 | N#N.O=C=CC=O                                                                                           | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -375.583  | -0.27844 | -0.07719 | 0.20125 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 22 | G16LOGParser | dsgdb9nsd_004015-3-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_004015-3-.log                                                      | .log          |       -1 |              1 | C[NH2+]C/C(C[O-])=N/[O-]                                                                               | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -417.832  | -0.0334  |  0.10345 | 0.13685 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 23 | G16LOGParser | dsgdb9nsd_004478-6-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_004478-6-.log                                                      | .log          |       -1 |              1 | Cc1[nH]c([O-])nc1O                                                                                     | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -415.566  | -0.00057 |  0.12524 | 0.12581 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 24 | G16LOGParser | dsgdb9nsd_004517-4-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_004517-4-.log                                                      | .log          |       -1 |              1 | Nc1[nH]c([O-])nc1O                                                                                     | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -431.606  |  0.00844 |  0.12415 | 0.11571 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 25 | G16LOGParser | dsgdb9nsd_004669-4.log                                                       | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_004669-4.log                                                       | .log          |        0 |              2 | O=CC1=C[N]C(O)=N1                                                                                      | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -414.256  | -0.28801 | -0.26252 | 0.02549 |      nan     |                  |      nan      |                   |   0.7589 |   0.5044 |
    | 26 | G16LOGParser | dsgdb9nsd_004738-2.log                                                       | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_004738-2.log                                                       | .log          |        0 |              2 | C#Cc1nc(O)[c]o1                                                                                        | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -396.863  | -0.26301 | -0.24718 | 0.01583 |      nan     |                  |      nan      |                   |   0.7794 |   0.5146 |
    | 27 | G16LOGParser | dsgdb9nsd_006075rearrange.log                                                | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_006075rearrange.log                                                | .log          |        0 |              1 | NC1(N)OCC(=O)O1                                                                                        | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -453.27   | -0.27969 | -0.01842 | 0.26127 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 28 | G16LOGParser | dsgdb9nsd_009986.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_009986.log                                                         | .log          |        0 |              1 | N#C[C-](C#N)/[NH+]=C/N                                                                                 | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -373.884  | -0.20306 | -0.07182 | 0.13124 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 29 | G16LOGParser | dsgdb9nsd_130366.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_130366.log                                                         | .log          |        0 |              1 | N[CH-]NC1=NC(=[OH+])N=N1                                                                               | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -466.411  | -0.19364 | -0.08163 | 0.11201 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 30 | G16LOGParser | dsgdb9nsd_131200-4-.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_131200-4-.log                                                      | .log          |       -1 |              1 | [O-]c1onc2c1CNC2                                                                                       | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -453.529  |  0.00581 |  0.07542 | 0.06961 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 31 | G16LOGParser | dsgdb9nsd_131200.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_131200.log                                                         | .log          |        0 |              1 | [O-]c1onc2c1C[NH2+]C2                                                                                  | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -454.13   | -0.18524 | -0.07291 | 0.11233 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 32 | G16LOGParser | dsgdb9nsd_131941-4+.log                                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_131941-4+.log                                                      | .log          |        1 |              1 | CN1[C+]=CC([N+](=O)[O-])=C1                                                                            | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -453.093  | -0.49069 | -0.38382 | 0.10687 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 33 | G16LOGParser | dsgdb9nsd_132072.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_132072.log                                                         | .log          |        0 |              1 | [NH2+]=c1[nH]ccc(F)c1[O-]                                                                              | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -478.256  | -0.196   | -0.05017 | 0.14583 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 34 | G16LOGParser | dsgdb9nsd_133826.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_133826.log                                                         | .log          |        0 |              1 | [C-]#C[C@@H]1N2C#[N+]C[C@]12C                                                                          | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -379.744  | -0.22785 | -0.1058  | 0.12205 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 35 | G16LOGParser | dsgdb9nsd_133858.log                                                         | /home/tmj/proj/MolOP/tests/test_files/g16log/dsgdb9nsd_133858.log                                                         | .log          |        0 |              1 | C[C@H]1C2C[C-]=C3C=[N+]1[C@@H]32                                                                       | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -364.836  | -0.18508 | -0.12276 | 0.06232 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 36 | G16LOGParser | molecule_0.log                                                               | /home/tmj/proj/MolOP/tests/test_files/g16log/molecule_0.log                                                               | .log          |        0 |              1 | C[C@@H]1[C@H](C)C2=C(CC2)[C@@H]1C                                                                      | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -390.536  | -0.21899 |  0.03946 | 0.25845 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 37 | G16LOGParser | r1_C2H3N3O_sp_g16.log                                                        | /home/tmj/proj/MolOP/tests/test_files/g16log/r1_C2H3N3O_sp_g16.log                                                        | .log          |        0 |              1 | C/C([O-])=N/[N+]#N                                                                                     | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -317.592  | -0.29408 | -0.06277 | 0.23131 |      nan     |                  |      nan      |                   | nan      | nan      |
    | 38 | G16LOGParser | r1_C8H9N3O2_sp_g16.log                                                       | /home/tmj/proj/MolOP/tests/test_files/g16log/r1_C8H9N3O2_sp_g16.log                                                       | .log          |        0 |              1 | CN/C([O-])=N/[O+]=N\c1ccccc1                                                                           | {'termination': 'Normal', 'SCF Done': True}                                                                                                                                                              | nan        | nan        | nan        | nan        |    nan    |   nan     |   nan     |   nan     |  -625.123  | -0.23807 | -0.12874 | 0.10933 |      nan     |                  |      nan      |                   | nan      | nan      |


我们关注这些日志文件中的 TS。


```python
TS_files = files.filter_TS()
print(TS_files.to_summary_df().to_markdown())
```

    |    | parser       | file_name                                                                    | file_path                                                                                                                 | file_format   |   charge |   multiplicity | SMILES                                                                                                 | status                                                                                                                                                                                                   |      ZPE |      TCE |      TCH |      TCG |   ZPE-Gas |    E-Gas |    H-Gas |    G-Gas |       sp |     HOMO |     LUMO |     GAP |   first freq | first freq tag   |   second freq | second freq tag   |   S**2 |   S |
    |---:|:-------------|:-----------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|:--------------|---------:|---------------:|:-------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------:|---------:|---------:|---------:|----------:|---------:|---------:|---------:|---------:|---------:|---------:|--------:|-------------:|:-----------------|--------------:|:------------------|-------:|----:|
    |  0 | G16LOGParser | RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log | /home/tmj/proj/MolOP/tests/test_files/g16log/RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log | .log          |       -1 |              1 | CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)CO2)[N-]1.C[CH-]C(=O)[N-]c1ccccc1.[Cu@TB7+5].[O-]c1ccccc1 | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': False, 'RMS     Displacement': True, 'termination': 'Normal', 'failure reason': '! Non-Optimized Parameters !', 'SCF Done': True} | 0.643965 | 0.682345 | 0.68329  | 0.572926 |  -1827.58 | -1827.54 | -1827.54 | -1827.65 | -1828.22 | -0.06449 |  0.0099  | 0.07439 |     -126.551 | True             |       14.3776 | False             |    nan | nan |
    |  1 | G16LOGParser | S_Ph_Ni_TS.log                                                               | /home/tmj/proj/MolOP/tests/test_files/g16log/S_Ph_Ni_TS.log                                                               | .log          |        0 |              1 | [Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1.c1c[nH]cn1.c1c[nH]cn1                                    | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': True, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                    | 0.361561 | 0.38868  | 0.389624 | 0.301123 |  -1548.69 | -1548.66 | -1548.66 | -1548.75 | -1549.05 | -0.18821 | -0.11033 | 0.07788 |     -185.866 | True             |       15.8954 | False             |     -0 |  -0 |
    |  2 | G16LOGParser | TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log                                      | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log                                      | .log          |        0 |              1 | CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(F)(F)F                                                   | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': False, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                   | 0.251882 | 0.277379 | 0.278323 | 0.192899 |  -1287.11 | -1287.08 | -1287.08 | -1287.17 | -1287.36 | -0.25863 | -0.10047 | 0.15816 |     -300.831 | True             |       16.957  | False             |    nan | nan |
    |  3 | G16LOGParser | TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log                                    | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log                                    | .log          |        0 |              1 | C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1                                                      | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': False, 'RMS     Displacement': True, 'termination': 'Normal', 'failure reason': '! Non-Optimized Parameters !', 'SCF Done': True} | 0.349175 | 0.373165 | 0.374109 | 0.294821 |  -6028.57 | -6028.55 | -6028.55 | -6028.63 | -6028.92 | -0.20404 | -0.04169 | 0.16235 |      -56.863 | True             |       29.9487 | False             |    nan | nan |
    |  4 | G16LOGParser | TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log                     | /home/tmj/proj/MolOP/tests/test_files/g16log/TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log                     | .log          |        0 |              1 | CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(F)(F)F                                                   | {'Maximum Force': True, 'RMS     Force': True, 'Maximum Displacement': True, 'RMS     Displacement': True, 'termination': 'Normal', 'SCF Done': True}                                                    | 0.252491 | 0.277515 | 0.27846  | 0.195753 |  -1287.12 | -1287.09 | -1287.09 | -1287.17 | -1287.37 | -0.25447 | -0.11005 | 0.14442 |     -252.41  | True             |       25.8663 | False             |    nan | nan |


原子数不多的情况


```python
TS_files[2][-1].rdmol
```




    
![png](read_file_and_extract_info_files/read_file_and_extract_info_9_0.png)
    



The numeric features (like totoal energy)
数值特征（如总能量）


```python
TS_files[2][-1].energy
```




-1287.35981 hartree/particle




```python
TS_files[2][-1].alpha_energy
```




    {'gap': 0.15816 <Unit('hartree / particle')>,
     'homo': -0.25863 <Unit('hartree / particle')>,
     'lumo': -0.10047 <Unit('hartree / particle')>}



输出单位可以转化为由[pint](https://pint.readthedocs.io/en/stable/)提供的其他单位。


```python
TS_files[2][-1].energy.to("kcal/mol")
```




-807830.4773030282 kilocalorie/mole



序列特征（如轨道、频率）


```python
TS_files[2][-1].alpha_FMO_orbits[:10]
```




    [-24.73197 <Unit('hartree / particle')>,
     -24.73103 <Unit('hartree / particle')>,
     -24.72836 <Unit('hartree / particle')>,
     -19.19597 <Unit('hartree / particle')>,
     -19.14519 <Unit('hartree / particle')>,
     -19.13671 <Unit('hartree / particle')>,
     -19.13657 <Unit('hartree / particle')>,
     -19.13465 <Unit('hartree / particle')>,
     -14.42477 <Unit('hartree / particle')>,
     -14.38013 <Unit('hartree / particle')>]



也提供无量纲特征


```python
TS_files[2][-1].dimensionless_alpha_energy
```




    {'gap': 0.15816, 'homo': -0.25863, 'lumo': -0.10047}



顺序无量纲特征作为生成器提供


```python
TS_files[2][-1].dimensionless_frequencies.__next__()
```




    {'freq': -300.8307,
     'is imaginary': True,
     'reduced masses': 11.3735,
     'IR intensities': 63.598,
     'force constants': 0.6064,
     'normal coordinates': [(0.08, -0.0, -0.05),
      (0.08, -0.05, -0.04),
      (0.06, -0.14, -0.0),
      (0.03, -0.04, -0.0),
      (0.03, -0.49, -0.01),
      (-0.04, 0.16, 0.03),
      (-0.13, -0.27, 0.03),
      (-0.03, -0.08, 0.03),
      (0.0, 0.01, 0.02),
      (-0.04, -0.01, -0.02),
      (-0.03, -0.07, -0.01),
      (-0.02, -0.02, 0.01),
      (-0.0, 0.03, 0.0),
      (0.08, 0.37, -0.02),
      (-0.08, 0.07, 0.0),
      (-0.04, -0.01, 0.0),
      (-0.12, 0.08, 0.01),
      (-0.04, 0.0, 0.0),
      (0.01, 0.54, 0.02),
      (0.06, 0.01, -0.0),
      (0.05, -0.05, -0.01),
      (0.06, -0.03, 0.01),
      (0.05, -0.02, 0.01),
      (0.09, 0.01, -0.05),
      (0.06, 0.01, -0.04),
      (0.1, -0.0, -0.05),
      (0.05, -0.17, 0.01),
      (0.03, 0.03, -0.03),
      (0.01, 0.03, 0.01),
      (-0.02, 0.04, 0.06),
      (0.0, 0.01, 0.01),
      (0.0, -0.02, -0.0),
      (-0.06, -0.0, 0.05),
      (0.07, -0.02, 0.01),
      (0.05, -0.01, 0.0),
      (0.04, -0.02, -0.0)]}



MolOP 可以从具有独特虚频方向的 TS 结构中推断出 TS 前和 TS 后结构。


```python
Draw.MolsToGridImage(TS_files[2][-1].possible_pre_post_ts())
```




    
![png](read_file_and_extract_info_files/read_file_and_extract_info_22_0.png)
    




```python
TS_files[2][-1].possible_pre_post_ts()[1]
```




    
![png](read_file_and_extract_info_files/read_file_and_extract_info_23_0.png)
    




```python
TS_files[2][-1].imaginary_frequencies
```




    [{'is imaginary': True,
      'freq': -300.8307 <Unit('reciprocal_centimeter')>,
      'reduced masses': 11.3735 <Unit('unified_atomic_mass_unit')>,
      'force constants': 0.6064 <Unit('millidyne / angstrom')>,
      'IR intensities': 63.598 <Unit('kilomole / mole')>,
      'normal coordinates': array([[ 0.08, -0.  , -0.05],
             [ 0.08, -0.05, -0.04],
             [ 0.06, -0.14, -0.  ],
             [ 0.03, -0.04, -0.  ],
             [ 0.03, -0.49, -0.01],
             [-0.04,  0.16,  0.03],
             [-0.13, -0.27,  0.03],
             [-0.03, -0.08,  0.03],
             [ 0.  ,  0.01,  0.02],
             [-0.04, -0.01, -0.02],
             [-0.03, -0.07, -0.01],
             [-0.02, -0.02,  0.01],
             [-0.  ,  0.03,  0.  ],
             [ 0.08,  0.37, -0.02],
             [-0.08,  0.07,  0.  ],
             [-0.04, -0.01,  0.  ],
             [-0.12,  0.08,  0.01],
             [-0.04,  0.  ,  0.  ],
             [ 0.01,  0.54,  0.02],
             [ 0.06,  0.01, -0.  ],
             [ 0.05, -0.05, -0.01],
             [ 0.06, -0.03,  0.01],
             [ 0.05, -0.02,  0.01],
             [ 0.09,  0.01, -0.05],
             [ 0.06,  0.01, -0.04],
             [ 0.1 , -0.  , -0.05],
             [ 0.05, -0.17,  0.01],
             [ 0.03,  0.03, -0.03],
             [ 0.01,  0.03,  0.01],
             [-0.02,  0.04,  0.06],
             [ 0.  ,  0.01,  0.01],
             [ 0.  , -0.02, -0.  ],
             [-0.06, -0.  ,  0.05],
             [ 0.07, -0.02,  0.01],
             [ 0.05, -0.01,  0.  ],
             [ 0.04, -0.02, -0.  ]], dtype=float16) <Unit('angstrom')>}]


