# Read files and extract information

There is an example of how to read files and extract information from them.

suppose you have some g16log files in a folder like this:


```python
! ls -l ../tests/test_files/g16log/
```

    total 33748
    -rw-r--r-- 1 tmj tmj  727184 Nov  4 20:14 11_Opt.log
    -rw-r--r-- 1 tmj tmj   53139 Jan  9 20:53 3_Sp.log
    -rw-r--r-- 1 tmj tmj 8008984 Feb 18 20:12 RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-Anion_From-IP_C-O-190_TS_Opt.log
    -rw-r--r-- 1 tmj tmj 1424205 Nov  4 20:14 S_Ph_Ni_TS.log
    -rw-r--r-- 1 tmj tmj 1910081 Jan 16 10:24 TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log
    -rw-r--r-- 1 tmj tmj   65746 Jan 16 10:28 TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log
    -rw-r--r-- 1 tmj tmj 3896697 Jan 15 22:37 TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log
    -rw-r--r-- 1 tmj tmj 1307536 Jan 16 10:28 TS_ts_guess_FaFxyx_template_4-18_6-13_optts_conf_g16.log
    -rw-r--r-- 1 tmj tmj  165117 Feb 17 22:09 dsgdb9nsd_000001-3-.log
    -rw-r--r-- 1 tmj tmj  823921 Jan 16 16:32 dsgdb9nsd_000107-3-.log
    -rw-r--r-- 1 tmj tmj  703600 Jan 16 14:55 dsgdb9nsd_000180-9-.log
    -rw-r--r-- 1 tmj tmj  385852 Jan 23 21:10 dsgdb9nsd_000484-1+.log
    -rw-r--r-- 1 tmj tmj  509600 Jan 22 21:03 dsgdb9nsd_000672-3+.log
    -rw-r--r-- 1 tmj tmj  583565 Jan 17 22:25 dsgdb9nsd_000696-4.log
    -rw-r--r-- 1 tmj tmj  380253 Jan 16 15:29 dsgdb9nsd_000763-2-.log
    -rw-r--r-- 1 tmj tmj  446864 Jan 22 18:22 dsgdb9nsd_000923-3+.log
    -rw-r--r-- 1 tmj tmj  523877 Feb 18 19:06 dsgdb9nsd_000955-3.log
    -rw-r--r-- 1 tmj tmj  383481 Jan 16 15:50 dsgdb9nsd_000958-3-.log
    -rw-r--r-- 1 tmj tmj  639463 Jan 16 15:07 dsgdb9nsd_001232-4-.log
    -rw-r--r-- 1 tmj tmj  535511 Jan 16 15:46 dsgdb9nsd_002924-8-.log
    -rw-r--r-- 1 tmj tmj  899961 Jan 18 10:54 dsgdb9nsd_003051-3.log
    -rw-r--r-- 1 tmj tmj  256920 Jan 10 20:09 dsgdb9nsd_003895.log
    -rw-r--r-- 1 tmj tmj  868095 Jan 16 18:07 dsgdb9nsd_004015-3-.log
    -rw-r--r-- 1 tmj tmj  442477 Jan 16 16:13 dsgdb9nsd_004478-6-.log
    -rw-r--r-- 1 tmj tmj  958696 Jan 16 16:06 dsgdb9nsd_004517-4-.log
    -rw-r--r-- 1 tmj tmj  672242 Jan 17 20:29 dsgdb9nsd_004669-4.log
    -rw-r--r-- 1 tmj tmj  453094 Jan 16 19:58 dsgdb9nsd_004738-2.log
    -rw-r--r-- 1 tmj tmj 1334751 Jan 10 20:23 dsgdb9nsd_006075rearrange.log
    -rw-r--r-- 1 tmj tmj  505226 Jan 10 22:06 dsgdb9nsd_009986.log
    -rw-r--r-- 1 tmj tmj  452946 Jan 10 21:53 dsgdb9nsd_130366.log
    -rw-r--r-- 1 tmj tmj  641086 Jan 17 20:29 dsgdb9nsd_131200-4-.log
    -rw-r--r-- 1 tmj tmj  499704 Jan 10 21:36 dsgdb9nsd_131200.log
    -rw-r--r-- 1 tmj tmj  911981 Jan 28 14:06 dsgdb9nsd_131941-4+.log
    -rw-r--r-- 1 tmj tmj  540614 Jan 10 21:45 dsgdb9nsd_132072.log
    -rw-r--r-- 1 tmj tmj  553016 Jan 10 21:43 dsgdb9nsd_133826.log
    -rw-r--r-- 1 tmj tmj  650970 Jan 10 20:54 dsgdb9nsd_133858.log
    -rw-r--r-- 1 tmj tmj  301634 Feb  1 18:21 molecule_0.log
    -rw-r--r-- 1 tmj tmj   22094 Jan 17 19:10 r1_C2H3N3O_sp_g16.log
    -rw-r--r-- 1 tmj tmj   38369 Jan 17 16:48 r1_C8H9N3O2_sp_g16.log


Looks like chaos? Yes, they are mixture from some different projects. MolOP can help you to read files and extract information from them through universal methods.


```python
from molop import AutoParser
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole  # for better drawing

IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 400, 400

files = AutoParser("../tests/test_files/g16log/*.log")
```

    MolOP parsing with 28 jobs: 100%|██████████| 39/39 [00:03<00:00, 10.88it/s]
    0 files failed to parse, 39 successfully parsed


We can first get the summary of the files. The summary contains the structure information(SMILES) and some key QM information of the files.


```python
files.to_summary_df()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>parser</th>
      <th>file_name</th>
      <th>file_path</th>
      <th>file_format</th>
      <th>charge</th>
      <th>multiplicity</th>
      <th>SMILES</th>
      <th>status</th>
      <th>ZPE</th>
      <th>TCE</th>
      <th>...</th>
      <th>sp</th>
      <th>HOMO</th>
      <th>LUMO</th>
      <th>GAP</th>
      <th>first freq</th>
      <th>first freq tag</th>
      <th>second freq</th>
      <th>second freq tag</th>
      <th>S**2</th>
      <th>S</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>G16LOGParser</td>
      <td>11_Opt.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/1...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>C=C[C@H]1[C@@H]2[C@H](C[C@@]1(C)O)OC[C@@H]2C</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.275585</td>
      <td>0.287903</td>
      <td>...</td>
      <td>-579.925317</td>
      <td>-0.21880</td>
      <td>-0.01452</td>
      <td>0.20428</td>
      <td>76.8500</td>
      <td>False</td>
      <td>94.3685</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>G16LOGParser</td>
      <td>3_Sp.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/3...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CCC[C@H]1CO[C@H]2C[C@@](C)(O)C[C@@H]12</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-581.780923</td>
      <td>-0.21593</td>
      <td>0.02581</td>
      <td>0.24174</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>G16LOGParser</td>
      <td>RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-An...</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/R...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)C...</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.643965</td>
      <td>0.682345</td>
      <td>...</td>
      <td>-1828.223960</td>
      <td>-0.06449</td>
      <td>0.00990</td>
      <td>0.07439</td>
      <td>-126.5511</td>
      <td>True</td>
      <td>14.3776</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>G16LOGParser</td>
      <td>S_Ph_Ni_TS.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/S...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>[Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1....</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.361561</td>
      <td>0.388680</td>
      <td>...</td>
      <td>-1549.053550</td>
      <td>-0.18821</td>
      <td>-0.11033</td>
      <td>0.07788</td>
      <td>-185.8656</td>
      <td>True</td>
      <td>15.8954</td>
      <td>False</td>
      <td>-0.0000</td>
      <td>-0.0000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>G16LOGParser</td>
      <td>TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(...</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.251882</td>
      <td>0.277379</td>
      <td>...</td>
      <td>-1287.359810</td>
      <td>-0.25863</td>
      <td>-0.10047</td>
      <td>0.15816</td>
      <td>-300.8307</td>
      <td>True</td>
      <td>16.9570</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>5</th>
      <td>G16LOGParser</td>
      <td>TS_4cGKps_ll_ad_4-18_6-13_sp_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(...</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-1288.848920</td>
      <td>-0.26367</td>
      <td>-0.10703</td>
      <td>0.15664</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>6</th>
      <td>G16LOGParser</td>
      <td>TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.349175</td>
      <td>0.373165</td>
      <td>...</td>
      <td>-6028.920709</td>
      <td>-0.20404</td>
      <td>-0.04169</td>
      <td>0.16235</td>
      <td>-56.8630</td>
      <td>True</td>
      <td>29.9487</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>7</th>
      <td>G16LOGParser</td>
      <td>TS_ts_guess_FaFxyx_template_4-18_6-13_optts_co...</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(...</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.252491</td>
      <td>0.277515</td>
      <td>...</td>
      <td>-1287.367861</td>
      <td>-0.25447</td>
      <td>-0.11005</td>
      <td>0.14442</td>
      <td>-252.4100</td>
      <td>True</td>
      <td>25.8663</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000001-3-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>[CH3-]</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-39.857502</td>
      <td>0.06578</td>
      <td>0.17809</td>
      <td>0.11231</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>9</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000107-3-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>[C-]#CCC#C</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-192.274958</td>
      <td>0.00382</td>
      <td>0.13824</td>
      <td>0.13442</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>10</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000180-9-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>CC/C(C)=N/[O-]</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-287.285392</td>
      <td>0.02146</td>
      <td>0.11362</td>
      <td>0.09216</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>11</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000484-1+.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>1</td>
      <td>1</td>
      <td>[C+]#CC#CC#C</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-228.646980</td>
      <td>-0.47931</td>
      <td>-0.43653</td>
      <td>0.04278</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>12</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000672-3+.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>1</td>
      <td>1</td>
      <td>[C+]1=C[C@H]2O[C@H]2C1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-268.397545</td>
      <td>-0.51954</td>
      <td>-0.32947</td>
      <td>0.19007</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>13</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000696-4.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>2</td>
      <td>C1[C@@H]2[C]3CN2[C@@H]13</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-248.755568</td>
      <td>-0.32784</td>
      <td>-0.18259</td>
      <td>0.14525</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>0.7535</td>
      <td>0.5018</td>
    </tr>
    <tr>
      <th>14</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000763-2-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>[H]/N=c1/[cH-]onn1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-316.936953</td>
      <td>0.00930</td>
      <td>0.13147</td>
      <td>0.12217</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>15</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000923-3+.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>1</td>
      <td>1</td>
      <td>O=C1C=C=NC=[NH+]1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-338.692192</td>
      <td>-0.49280</td>
      <td>-0.35012</td>
      <td>0.14268</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>16</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000955-3.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>2</td>
      <td>Oc1nc[c]cn1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-338.979788</td>
      <td>-0.28326</td>
      <td>-0.27967</td>
      <td>0.00359</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>0.7570</td>
      <td>0.5035</td>
    </tr>
    <tr>
      <th>17</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_000958-3-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>Oc1n[c-]ncn1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-355.093513</td>
      <td>-0.00299</td>
      <td>0.12011</td>
      <td>0.12310</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>18</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_001232-4-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>N/C=N/C=C(/[O-])O</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-377.444816</td>
      <td>0.00144</td>
      <td>0.10845</td>
      <td>0.10701</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>19</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_002924-8-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>COC1=NCC[N-]1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-341.474822</td>
      <td>0.01000</td>
      <td>0.12531</td>
      <td>0.11531</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>20</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_003051-3.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>2</td>
      <td>O=CN[CH][C@H]1CO1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-361.272657</td>
      <td>-0.28613</td>
      <td>-0.19199</td>
      <td>0.09414</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>0.7565</td>
      <td>0.5032</td>
    </tr>
    <tr>
      <th>21</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_003895.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>N#N.O=C=CC=O</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-375.582679</td>
      <td>-0.27844</td>
      <td>-0.07719</td>
      <td>0.20125</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>22</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_004015-3-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>C[NH2+]C/C(C[O-])=N/[O-]</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-417.832156</td>
      <td>-0.03340</td>
      <td>0.10345</td>
      <td>0.13685</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>23</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_004478-6-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>Cc1[nH]c([O-])nc1O</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-415.565992</td>
      <td>-0.00057</td>
      <td>0.12524</td>
      <td>0.12581</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>24</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_004517-4-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>Nc1[nH]c([O-])nc1O</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-431.605850</td>
      <td>0.00844</td>
      <td>0.12415</td>
      <td>0.11571</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>25</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_004669-4.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>2</td>
      <td>O=CC1=C[N]C(O)=N1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-414.256471</td>
      <td>-0.28801</td>
      <td>-0.26252</td>
      <td>0.02549</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>0.7589</td>
      <td>0.5044</td>
    </tr>
    <tr>
      <th>26</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_004738-2.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>2</td>
      <td>C#Cc1nc(O)[c]o1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-396.862586</td>
      <td>-0.26301</td>
      <td>-0.24718</td>
      <td>0.01583</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>0.7794</td>
      <td>0.5146</td>
    </tr>
    <tr>
      <th>27</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_006075rearrange.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>NC1(N)OCC(=O)O1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-453.270160</td>
      <td>-0.27969</td>
      <td>-0.01842</td>
      <td>0.26127</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>28</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_009986.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>N#C[C-](C#N)/[NH+]=C/N</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-373.883900</td>
      <td>-0.20306</td>
      <td>-0.07182</td>
      <td>0.13124</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>29</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_130366.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>N[CH-]NC1=NC(=[OH+])N=N1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-466.410598</td>
      <td>-0.19364</td>
      <td>-0.08163</td>
      <td>0.11201</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>30</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_131200-4-.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>[O-]c1onc2c1CNC2</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-453.529216</td>
      <td>0.00581</td>
      <td>0.07542</td>
      <td>0.06961</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>31</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_131200.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>[O-]c1onc2c1C[NH2+]C2</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-454.129865</td>
      <td>-0.18524</td>
      <td>-0.07291</td>
      <td>0.11233</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>32</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_131941-4+.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>1</td>
      <td>1</td>
      <td>CN1[C+]=CC([N+](=O)[O-])=C1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-453.092643</td>
      <td>-0.49069</td>
      <td>-0.38382</td>
      <td>0.10687</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>33</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_132072.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>[NH2+]=c1[nH]ccc(F)c1[O-]</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-478.256385</td>
      <td>-0.19600</td>
      <td>-0.05017</td>
      <td>0.14583</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>34</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_133826.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>[C-]#C[C@@H]1N2C#[N+]C[C@]12C</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-379.743730</td>
      <td>-0.22785</td>
      <td>-0.10580</td>
      <td>0.12205</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>35</th>
      <td>G16LOGParser</td>
      <td>dsgdb9nsd_133858.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/d...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>C[C@H]1C2C[C-]=C3C=[N+]1[C@@H]32</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-364.835805</td>
      <td>-0.18508</td>
      <td>-0.12276</td>
      <td>0.06232</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>36</th>
      <td>G16LOGParser</td>
      <td>molecule_0.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/m...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>C[C@@H]1[C@H](C)C2=C(CC2)[C@@H]1C</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-390.536199</td>
      <td>-0.21899</td>
      <td>0.03946</td>
      <td>0.25845</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>37</th>
      <td>G16LOGParser</td>
      <td>r1_C2H3N3O_sp_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/r...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>C/C([O-])=N/[N+]#N</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-317.592367</td>
      <td>-0.29408</td>
      <td>-0.06277</td>
      <td>0.23131</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>38</th>
      <td>G16LOGParser</td>
      <td>r1_C8H9N3O2_sp_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/r...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CN/C([O-])=N/[O+]=N\c1ccccc1</td>
      <td>{'termination': 'Normal', 'SCF Done': True}</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>-625.122643</td>
      <td>-0.23807</td>
      <td>-0.12874</td>
      <td>0.10933</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>39 rows × 26 columns</p>
</div>



We concern about the TS in those log files.


```python
TS_files = files.filter_TS()
TS_files.to_summary_df()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>parser</th>
      <th>file_name</th>
      <th>file_path</th>
      <th>file_format</th>
      <th>charge</th>
      <th>multiplicity</th>
      <th>SMILES</th>
      <th>status</th>
      <th>ZPE</th>
      <th>TCE</th>
      <th>...</th>
      <th>sp</th>
      <th>HOMO</th>
      <th>LUMO</th>
      <th>GAP</th>
      <th>first freq</th>
      <th>first freq tag</th>
      <th>second freq</th>
      <th>second freq tag</th>
      <th>S**2</th>
      <th>S</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>G16LOGParser</td>
      <td>RE_BOX-Anion-Real_Cu-III-Phenol_Major-Amide-An...</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/R...</td>
      <td>.log</td>
      <td>-1</td>
      <td>1</td>
      <td>CC(C)(C)[C@@H]1COC(=[C-]C2[N-][C@H](C(C)(C)C)C...</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.643965</td>
      <td>0.682345</td>
      <td>...</td>
      <td>-1828.223960</td>
      <td>-0.06449</td>
      <td>0.00990</td>
      <td>0.07439</td>
      <td>-126.5511</td>
      <td>True</td>
      <td>14.3776</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>G16LOGParser</td>
      <td>S_Ph_Ni_TS.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/S...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>[Br-].[Br-].[Ni@OH26+3]c1ccccc1.[S-]Cc1ccccc1....</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.361561</td>
      <td>0.388680</td>
      <td>...</td>
      <td>-1549.053550</td>
      <td>-0.18821</td>
      <td>-0.11033</td>
      <td>0.07788</td>
      <td>-185.8656</td>
      <td>True</td>
      <td>15.8954</td>
      <td>False</td>
      <td>-0.0</td>
      <td>-0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>G16LOGParser</td>
      <td>TS_4cGKps_ll_ad_4-18_6-13_optts_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CNC(=O)C#[N+]/C(C(C)=O)=C(/C)[O-].COC(=O)C#CC(...</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.251882</td>
      <td>0.277379</td>
      <td>...</td>
      <td>-1287.359810</td>
      <td>-0.25863</td>
      <td>-0.10047</td>
      <td>0.15816</td>
      <td>-300.8307</td>
      <td>True</td>
      <td>16.9570</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>G16LOGParser</td>
      <td>TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.log</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>C=[N+](C)[N-]C.COC(=O)[C@@]1(OC)C#CC(Br)(Br)CCCC1</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.349175</td>
      <td>0.373165</td>
      <td>...</td>
      <td>-6028.920709</td>
      <td>-0.20404</td>
      <td>-0.04169</td>
      <td>0.16235</td>
      <td>-56.8630</td>
      <td>True</td>
      <td>29.9487</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>G16LOGParser</td>
      <td>TS_ts_guess_FaFxyx_template_4-18_6-13_optts_co...</td>
      <td>/home/tmj/proj/MolOP/tests/test_files/g16log/T...</td>
      <td>.log</td>
      <td>0</td>
      <td>1</td>
      <td>CNC(=O)C#[N+]/C(C(C)=O)=C(\C)[O-].COC(=O)C#CC(...</td>
      <td>{'Maximum Force': True, 'RMS     Force': True,...</td>
      <td>0.252491</td>
      <td>0.277515</td>
      <td>...</td>
      <td>-1287.367861</td>
      <td>-0.25447</td>
      <td>-0.11005</td>
      <td>0.14442</td>
      <td>-252.4100</td>
      <td>True</td>
      <td>25.8663</td>
      <td>False</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 26 columns</p>
</div>



A case with not too many atoms


```python
TS_files[2][-1].rdmol
```




    
![png](read_file_and_extract_info_files/read_file_and_extract_info_9_0.png)
    



The numeric features (like totoal energy)


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



The unit can be transformed to other units powered by [pint](https://pint.readthedocs.io/en/stable/).


```python
TS_files[2][-1].energy.to("kcal/mol")
```




-807830.4773030282 kilocalorie/mole



The sequential features (like orbitals, freqs)


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



Also offer dimensionless features


```python
TS_files[2][-1].dimensionless_alpha_energy
```




    {'gap': 0.15816, 'homo': -0.25863, 'lumo': -0.10047}



The sequential dimensionless features are provided as generator


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



MolOP can infer the pre- and post- TS structures from the TS structure with its unique imaginary frequency orientation.


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


