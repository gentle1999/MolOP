# Structure recovery cases


```python
from molop.io import AutoParser
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 400, 400
```

## molecule_0.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/molecule_0.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 653.11it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_3_1.png)
    



## dsgdb9nsd_131941-4+.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_131941-4+.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 83.35it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_5_1.png)
    




```python
test_file[0][-1].rdmol_no_conformer
```




    
![png](structure_recovery_cases_files/structure_recovery_cases_6_0.png)
    



## RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 34.55it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_8_1.png)
    



## add_1.fchk


```python
test_file = AutoParser(
    "../../tests/test_files/g16fchk/add_1.fchk",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  5.03it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_10_1.png)
    



## add_0.fchk


```python
test_file = AutoParser(
    "../../tests/test_files/g16fchk/add_0.fchk",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 73.23it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_12_1.png)
    



## dsgdb9nsd_000007-6.fchk


```python
test_file = AutoParser(
    "../../tests/test_files/g16fchk/dsgdb9nsd_000007-6.fchk",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 285.81it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_14_1.png)
    



## dsgdb9nsd_000001-3.fchk


```python
test_file = AutoParser(
    "../../tests/test_files/g16fchk/dsgdb9nsd_000001-3.fchk",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 634.54it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_16_1.png)
    



## dsgdb9nsd_000484-1+.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000484-1+.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 113.33it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_18_1.png)
    



## dsgdb9nsd_000672-3+.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000672-3+.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 138.28it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_20_1.png)
    



## dsgdb9nsd_000923-3+.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000923-3+.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 130.81it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_22_1.png)
    



## dsgdb9nsd_003051-3.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_003051-3.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 69.03it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_24_1.png)
    



## dsgdb9nsd_000696-4.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000696-4.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 72.59it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_26_1.png)
    



## dsgdb9nsd_004669-4.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_004669-4.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 56.04it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_28_1.png)
    



## dsgdb9nsd_131200-4-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_131200-4-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 91.51it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_30_1.png)
    



## dsgdb9nsd_000955-3.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000955-3.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 67.79it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_32_1.png)
    




```python
test_file[0][-1].rdmol_no_conformer
```




    
![png](structure_recovery_cases_files/structure_recovery_cases_33_0.png)
    



## r1_C2H3N3O_sp_g16.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/r1_C2H3N3O_sp_g16.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 1503.33it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_35_1.png)
    



## r1_C8H9N3O2_sp_g16.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/r1_C8H9N3O2_sp_g16.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 1095.40it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_37_1.png)
    



## dsgdb9nsd_004738-2.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_004738-2.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 64.48it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_39_1.png)
    



## dsgdb9nsd_004015-3-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_004015-3-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 104.88it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_41_1.png)
    



## dsgdb9nsd_000107-3-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000107-3-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 161.34it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_43_1.png)
    



## dsgdb9nsd_004517-4-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_004478-6-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 107.92it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_45_1.png)
    




```python
test_file[0][-1].rdmol_no_conformer
```




    
![png](structure_recovery_cases_files/structure_recovery_cases_46_0.png)
    



## dsgdb9nsd_004517-4-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_004517-4-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 75.57it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_48_1.png)
    



## dsgdb9nsd_000958-3-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000958-3-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 113.22it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_50_1.png)
    




```python
test_file[0][-1].rdmol_no_conformer
```




    
![png](structure_recovery_cases_files/structure_recovery_cases_51_0.png)
    



## dsgdb9nsd_002924-8-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_002924-8-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 97.61it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_53_1.png)
    



## dsgdb9nsd_000763-2-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000763-2-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 158.35it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_55_1.png)
    



## dsgdb9nsd_001232-4-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_001232-4-.log",
    n_jobs=1,
    only_extract_structure=True,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 127.49it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_57_1.png)
    



## dsgdb9nsd_000180-9-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000180-9-.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  9.63it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_59_1.png)
    



## TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 50.25it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_61_1.png)
    



## TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log",
    n_jobs=1,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 20.39it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_63_1.png)
    



## TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log",
    n_jobs=1,
    only_last_frame=True,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 20.97it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_65_1.png)
    



## TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  2.33it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_67_1.png)
    



## irc.out


```python
test_file = AutoParser(
    "../../tests/test_files/g16irc/irc.out",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  7.19it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_69_1.png)
    



## 3_Sp.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/3_Sp.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 70.30it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_71_1.png)
    



## 11_Opt.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/11_Opt.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  8.18it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_73_1.png)
    



## dsgdb9nsd_003895.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_003895.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 27.73it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_75_1.png)
    



## dsgdb9nsd_006075rearrange.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_006075rearrange.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  8.22it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_77_1.png)
    



## dsgdb9nsd_130366.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_130366.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 16.99it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_79_1.png)
    



## dsgdb9nsd_131200.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_131200.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 15.64it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_81_1.png)
    



## dsgdb9nsd_132072.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_132072.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 14.90it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_83_1.png)
    



## dsgdb9nsd_133826.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_133826.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 15.18it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_85_1.png)
    



## dsgdb9nsd_133858.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_133858.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 13.07it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_87_1.png)
    



## S_Ph_Ni_TS.gjf


```python
test_file = AutoParser(
    "../../tests/test_files/g16gjf/S_Ph_Ni_TS.gjf",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 2091.92it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_89_1.png)
    



## S_Ph_Ni_TS.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/S_Ph_Ni_TS.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00,  5.81it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_91_1.png)
    



## test.gjf


```python
test_file = AutoParser(
    "../../tests/test_files/g16gjf/test.gjf",
    n_jobs=1,
    charge=0,
    multiplicity=3,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 1633.93it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_93_1.png)
    



## test.xyz


```python
test_file = AutoParser(
    "../../tests/test_files/xyz/test.xyz",
    n_jobs=1,
    charge=0,
    multiplicity=3,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 1249.79it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_95_1.png)
    



## dsgdb9nsd_009986.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_009986.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 16.66it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_97_1.png)
    



## dsgdb9nsd_000001-3-.log


```python
test_file = AutoParser(
    "../../tests/test_files/g16log/dsgdb9nsd_000001-3-.log",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 60.34it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_99_1.png)
    



## xtb_6_3_2_opt.out
Out file is the printout on screen of the xtb calculation. 
you can get it by:

```bash
xtb ***.xyz > ***.out
```


```python
test_file = AutoParser(
    "../../tests/test_files/xtbout/xtb_6_3_2_opt.out",
    n_jobs=1,
)
test_file[0][-1].rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 202.38it/s]
    0 files failed to parse, 1 successfully parsed





    
![png](structure_recovery_cases_files/structure_recovery_cases_101_1.png)
    


