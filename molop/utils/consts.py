k = 8.9875517923e9
qe = 1.602176634 * 10**-19


"""
The cutoff_epsilon_sigma_dict is a copy of [LennardJones612_UniversalShifted.params](https://openkim.org/files/MO_959249795837_003/LennardJones612_UniversalShifted.params)
"""
cutoff_epsilon_sigma_dict = {
    "H": (2.2094300, 4.4778900, 0.5523570),
    "He": (1.9956100, 0.0009421, 0.4989030),
    "Li": (9.1228000, 1.0496900, 2.2807000),
    "Be": (6.8421000, 0.5729420, 1.7105300),
    "B": (6.0581100, 2.9670300, 1.5145300),
    "C": (5.4166600, 6.3695300, 1.3541700),
    "N": (5.0603000, 9.7537900, 1.2650800),
    "O": (4.7039500, 5.1264700, 1.1759900),
    "F": (4.0625000, 1.6059200, 1.0156200),
    "Ne": (4.1337700, 0.0036471, 1.0334400),
    "Na": (11.8311000, 0.7367450, 2.9577800),
    "Mg": (10.0493000, 0.0785788, 2.5123300),
    "Al": (8.6239000, 2.7006700, 2.1559700),
    "Si": (7.9111800, 3.1743100, 1.9778000),
    "P": (7.6260900, 5.0305000, 1.9065200),
    "S": (7.4835500, 4.3692700, 1.8708900),
    "Cl": (7.2697300, 4.4832800, 1.8174300),
    "Ar": (7.5548200, 0.0123529, 1.8887100),
    "K": (14.4682000, 0.5517990, 3.6170500),
    "Ca": (12.5439000, 0.1326790, 3.1359600),
    "Sc": (12.1162000, 1.6508000, 3.0290600),
    "Ti": (11.4035000, 1.1802700, 2.8508800),
    "V": (10.9046000, 2.7524900, 2.7261500),
    "Cr": (9.9067900, 1.5367900, 2.4767000),
    "Mn": (9.9067900, 0.5998880, 2.4767000),
    "Fe": (9.4078900, 1.1844200, 2.3519700),
    "Co": (8.9802600, 1.2776900, 2.2450600),
    "Ni": (8.8377200, 2.0757200, 2.2094300),
    "Cu": (9.4078900, 2.0446300, 2.3519700),
    "Zn": (8.6951700, 0.1915460, 2.1737900),
    "Ga": (8.6951700, 1.0642000, 2.1737900),
    "Ge": (8.5526300, 2.7017100, 2.1381600),
    "As": (8.4813600, 3.9599000, 2.1203400),
    "Se": (8.5526300, 3.3867700, 2.1381600),
    "Br": (8.5526300, 1.9706300, 2.1381600),
    "Kr": (8.2675400, 0.0173276, 2.0668900),
    "Rb": (15.6798000, 0.4682650, 3.9199500),
    "Sr": (13.8980000, 0.1339230, 3.4745100),
    "Y": (13.5417000, 2.7597500, 3.3854200),
    "Zr": (12.4726000, 3.0520100, 3.1181500),
    "Nb": (11.6886000, 5.2782000, 2.9221500),
    "Mo": (10.9759000, 4.4749900, 2.7439700),
    "Tc": (10.4770000, 3.3815900, 2.6192400),
    "Ru": (10.4057000, 1.9617200, 2.6014200),
    "Rh": (10.1206000, 2.4058200, 2.5301500),
    "Pd": (9.9067900, 1.3709700, 2.4767000),
    "Ag": (10.3344000, 1.6497600, 2.5836100),
    "Cd": (10.2632000, 0.0377447, 2.5657900),
    "In": (10.1206000, 0.8113140, 2.5301500),
    "Sn": (9.9067900, 1.9005700, 2.4767000),
    "Sb": (9.9067900, 3.0882800, 2.4767000),
    "Te": (9.8355200, 2.6312300, 2.4588800),
    "I": (9.9067900, 1.5393800, 2.4767000),
    "Xe": (9.9780700, 0.0238880, 2.4945200),
    "Cs": (17.3903000, 0.4166420, 4.3475900),
    "Ba": (15.3235000, 1.9000000, 3.8308600),
    "La": (14.7533000, 2.4996100, 3.6883200),
    "Ce": (14.5395000, 2.5700800, 3.6348700),
    "Pr": (14.4682000, 1.2994600, 3.6170500),
    "Nd": (14.3257000, 0.8196050, 3.5814100),
    "Pm": (14.1831000, 3.2413400, 3.5457800),
    "Sm": (14.1118000, 0.5211220, 3.5279600),
    "Eu": (14.1118000, 0.4299180, 3.5279600),
    "Gd": (13.9693000, 2.0995600, 3.4923200),
    "Tb": (13.8267000, 1.3999900, 3.4566900),
    "Dy": (13.6842000, 0.6900550, 3.4210500),
    "Ho": (13.6842000, 0.6900550, 3.4210500),
    "Er": (13.4704000, 0.7387660, 3.3676000),
    "Tm": (13.5417000, 0.5211220, 3.3854200),
    "Yb": (13.3278000, 0.1303990, 3.3319600),
    "Lu": (13.3278000, 1.4331500, 3.3319600),
    "Hf": (12.4726000, 3.3608600, 3.1181500),
    "Ta": (12.1162000, 4.0034300, 3.0290600),
    "W": (11.5460000, 6.8638900, 2.8865100),
    "Re": (10.7621000, 4.4387100, 2.6905100),
    "Os": (10.2632000, 4.2625300, 2.5657900),
    "Ir": (10.0493000, 3.7028700, 2.5123300),
    "Pt": (9.6929800, 3.1401000, 2.4232400),
    "Au": (9.6929800, 2.3058000, 2.4232400),
    "Hg": (9.4078900, 0.0454140, 2.3519700),
    "Tl": (10.3344000, 0.5770870, 2.5836100),
    "Pb": (10.4057000, 0.8589880, 2.6014200),
    "Bi": (10.5482000, 2.0798700, 2.6370600),
    "Po": (9.9780700, 1.8995300, 2.4945200),
    "At": (10.6908000, 1.3854420, 2.6727000),
    "Rn": (10.6908000, 0.0214992, 2.6727000),
    "Fr": (18.5307000, 0.3749778, 4.6326700),
    "Ra": (15.7511000, 1.7100000, 3.9377700),
    "Ac": (15.3235000, 2.2496490, 3.8308600),
    "Th": (14.6820000, 2.3130720, 3.6705000),
    "Pa": (14.2544000, 1.1695140, 3.5635900),
    "U": (13.9693000, 0.7376445, 3.4923200),
    "Np": (13.5417000, 2.9172060, 3.3854200),
    "Pu": (13.3278000, 0.4690098, 3.3319600),
    "Am": (12.8289000, 0.3869262, 3.2072400),
    "Cm": (12.0450000, 1.8896040, 3.0112400),
    "Bk": (11.9737000, 1.2599910, 2.9934200),
    "Cf": (11.9737000, 0.6210495, 2.9934200),
    "Es": (11.7599000, 0.6210495, 2.9399700),
    "Fm": (11.9024000, 0.6648894, 2.9756000),
    "Md": (12.3300000, 0.4690098, 3.0825100),
    "No": (12.5439000, 0.1173591, 3.1359600),
    "Lr": (11.4748000, 1.2898350, 2.8686900),
    "Rf": (11.1897000, 3.0247740, 2.7974200),
    "Db": (10.6195000, 3.6030870, 2.6548800),
    "Sg": (10.1919000, 6.1775010, 2.5479700),
    "Bh": (10.0493000, 3.9948390, 2.5123300),
    "Hs": (9.5504300, 3.8362770, 2.3876100),
    "Mt": (9.1940700, 3.3325830, 2.2985200),
    "Ds": (9.1228000, 2.8260900, 2.2807000),
    "Rg": (8.6239000, 2.0752200, 2.1559700),
    "Cn": (8.6951700, 0.0408726, 2.1737900),
    "Nh": (9.6929800, 0.5193783, 2.4232400),
    "Fl": (10.1919000, 0.7730892, 2.5479700),
    "Mc": (11.5460000, 1.8718830, 2.8865100),
    "Lv": (12.4726000, 1.7095770, 3.1181500),
    "Ts": (11.7599000, 1.2468978, 2.9399700),
    "Og": (11.1897000, 0.0193493, 2.7974200),
}
