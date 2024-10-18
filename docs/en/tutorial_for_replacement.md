# Tutorial for Replacement

Replacing substituents on the 3D structure of a molecule is a very common requirement. For small numbers of substitutions, manual modifications can be made using visualization tools such as GaussView. An automated tool becomes necessary when faced with large-scale dataset expansion.

RDKit provides functions for substituent replacement, which perform well on 2D molecular maps but are difficult to maintain robustness on 3D structural substitutions. A set of automatic replacement functions for single-bond linked substituents has been implemented in MolOP, allowing to keep the molecular backbone and the substituents conformationally invariant to each other and to find a suitable angle of rotation making the site resistance small.


## Demo Data

A C-H activation reaction of a ferrocene derivative catalyzed by an amino acid amide, where a possible transition state is shown in the figure. Note that the hydrogen atom on number 17 is activated.


```python
from molop import AutoParser
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 600, 600

template_rp = AutoParser("../../tests/test_files/g16gjf/Template-Rp.gjf")[0][0]
template_rp.rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 1104.93it/s]
    INFO - 0 files failed to parse, 1 successfully parsed





    
![png](tutorial_for_replacement_files/tutorial_for_replacement_2_1.png)
    



The activity of this reaction may be affected by the derivatization of ferrocene itself and the derivatization of amino acid amide. Therefore, we hope to modify this template to produce the corresponding transition state in batches.


```python
import pandas as pd
from tqdm import tqdm

tqdm.pandas()

demo_data = pd.read_csv("../data/replacement_demo_data.csv", index_col=0)
demo_data
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
      <th>Catalyst</th>
      <th>Reactant</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CC(C)(C)OC(=O)N[C@H](C(=O)O)C(C)(C)C</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>5</th>
      <td>CC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>6</th>
      <td>CC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>7</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>8</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>9</th>
      <td>CC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>10</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>11</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>12</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>13</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>14</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>15</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>16</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>17</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>18</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>19</th>
      <td>CC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>20</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>21</th>
      <td>CC(C)C[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>22</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7(Br)C28</td>
    </tr>
    <tr>
      <th>23</th>
      <td>CC(=O)N[C@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>24</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>25</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>26</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>27</th>
      <td>CC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>28</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7(P(=O)(c1...</td>
    </tr>
    <tr>
      <th>29</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7(C(O)(c1c...</td>
    </tr>
    <tr>
      <th>30</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>31</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>32</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>33</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>C1CCN(CC23C4C5C6C2[Fe]56432789C3C2C7C8C39)C1</td>
    </tr>
    <tr>
      <th>34</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>35</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>36</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>C1CCN(CC23C4C5C6C2[Fe]56432789C3C2C7C8C39)C1</td>
    </tr>
    <tr>
      <th>37</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>38</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>39</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>40</th>
      <td>CC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>41</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>42</th>
      <td>CC(C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>43</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>44</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>45</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>46</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>47</th>
      <td>CC[C@H](C)[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>48</th>
      <td>CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
    <tr>
      <th>49</th>
      <td>C[C@H](NC(=O)OC(C)(C)C)C(=O)O</td>
      <td>CN(C)CC12C3C4C5C1[Fe]45321678C2C1C6C7C28</td>
    </tr>
  </tbody>
</table>
</div>



Sometimes we can tell what substituent should be replaced at what position, but more often than not, we only know the structure of the substrate. The demo data does exactly that. We can quickly access the substituents to be replaced with the help of RDKit and known substrate common parts.


```python
from rdkit import Chem
from rdkit.Chem.rdmolops import ReplaceCore


def trasform_substituent(smiles):
    rwmol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    rwmol.RemoveAtom(0)
    rwmol.GetAtomWithIdx(0).SetNumRadicalElectrons(
        rwmol.GetAtomWithIdx(0).GetNumRadicalElectrons() + 1
    )
    return Chem.MolToSmiles(rwmol.GetMol())


def transform_substituents(row: pd.Series):
    smis = row[0].split(".")
    return pd.Series([trasform_substituent(s) for s in smis])


core_fe = Chem.MolFromSmiles("CC12C3C4C5C1[Fe]23451678C2C1C6C7C28")
core_animo = Chem.MolFromSmiles("O=CNCC(=O)O")

fe_substituents = demo_data.apply(
    lambda row: pd.Series(
        Chem.MolToSmiles(
            ReplaceCore(
                Chem.MolFromSmiles(row["Reactant"]),
                core_fe,
            )
        )
    ),
    axis=1,
)

amino_substituents = demo_data.apply(
    lambda row: pd.Series(
        Chem.MolToSmiles(
            ReplaceCore(
                Chem.MolFromSmiles(row["Catalyst"]),
                core_animo,
            )
        )
    ),
    axis=1,
)
```


```python
fe_substituents_df = fe_substituents.apply(transform_substituents, axis=1).fillna("[H]")
fe_substituents_df.columns = ["fe_substituents_0", "fe_substituents_1"]
amino_substituents_df = amino_substituents.apply(transform_substituents, axis=1).fillna(
    "[H]"
)
amino_substituents_df.columns = ["amino_substituents_0", "amino_substituents_1"]
replacemnet_df = pd.concat([fe_substituents_df, amino_substituents_df], axis=1)
replacemnet_df
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
      <th>fe_substituents_0</th>
      <th>fe_substituents_1</th>
      <th>amino_substituents_0</th>
      <th>amino_substituents_1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>3</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[C](C)C</td>
    </tr>
    <tr>
      <th>5</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>8</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>9</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>12</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>13</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>14</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>15</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>16</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>17</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>18</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>19</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>20</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>21</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]C(C)C</td>
    </tr>
    <tr>
      <th>22</th>
      <td>C[N]C</td>
      <td>[Br]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>23</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>24</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>25</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>26</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>27</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>28</th>
      <td>C[N]C</td>
      <td>O=[P](c1ccccc1)c1ccccc1</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>29</th>
      <td>C[N]C</td>
      <td>O[C](c1ccccc1)c1ccccc1</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>30</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>31</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>32</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>33</th>
      <td>C1CC[N]C1</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>34</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>35</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>36</th>
      <td>C1CC[N]C1</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>37</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>38</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>39</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>40</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>[CH3]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>41</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>42</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]C</td>
    </tr>
    <tr>
      <th>43</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>44</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>45</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>46</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>47</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>C[CH]CC</td>
    </tr>
    <tr>
      <th>48</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH2]c1ccccc1</td>
    </tr>
    <tr>
      <th>49</th>
      <td>C[N]C</td>
      <td>[H]</td>
      <td>CC(C)(C)[O]</td>
      <td>[CH3]</td>
    </tr>
  </tbody>
</table>
</div>



Subsequently, substitutions are made at the specified positions (based on the `bind_idx` parameter) according to the obtained substituents. Note that substituent substitution only guarantees the absolute atom number of the substructure before `bind_idx`; the atom number of the new substituent will come last. Therefore, it is necessary to ensure that positions with a larger `bind_idx` are substituted first.


```python
def preprocess_template_1(row: pd.Series):
    ts_like_rp = (
        template_rp.replace_substituent(
            "[H]", row["fe_substituents_1"], bind_idx=34, angle_split=20
        )
        .replace_substituent(
            "[H]", row["fe_substituents_0"], bind_idx=32, angle_split=20
        )
        .replace_substituent(
            "[H]", row["amino_substituents_1"], bind_idx=29, angle_split=20
        )
        .replace_substituent("[H]", row["amino_substituents_0"], 28, angle_split=20)
    )
    return ts_like_rp.rdmol


replacemnet_rdmols = replacemnet_df.progress_apply(preprocess_template_1, axis=1)
```

    100%|██████████| 50/50 [00:10<00:00,  4.81it/s]



```python
replacemnet_rdmols.iloc[0]
```




    
![png](tutorial_for_replacement_files/tutorial_for_replacement_10_0.png)
    



If the representation of the new substituent is based on SMILES, then its conformation is generated by RDKit, which means that its structural accuracy is not higher than the force field level. In addition, RDKit's conformational generation algorithm does not support metals very well, in particular it is difficult to generate the correct coordination conformer. In this scenario, you can just pass the rdmol object and just make sure that there is a free radical in it, or if there is more than one then you can specify it by the relative order within. 


```python
substituent = AutoParser("../../tests/test_files/xyz/dsgdb9nsd_046611-4/12.xyz")[0][
    0
].rdmol

template_rp.replace_substituent("[H]", substituent, bind_idx=34).rdmol
```

    MolOP parsing with single thread: 100%|██████████| 1/1 [00:00<00:00, 1461.94it/s]
    INFO - 0 files failed to parse, 1 successfully parsed





    
![png](tutorial_for_replacement_files/tutorial_for_replacement_12_1.png)
    




```python

```
