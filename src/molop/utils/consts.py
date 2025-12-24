from dataclasses import dataclass


K = 8.9875517923e9
QE = 1.602176634 * 10**-19


"""
The cutoff_epsilon_sigma_dict is a copy of
[LennardJones612_UniversalShifted.params](https://openkim.org/files/MO_959249795837_003/LennardJones612_UniversalShifted.params)
"""


@dataclass
class LennardJonesParams:
    epsilon: float
    sigma: float
    cutoff: float


CUTOFF_EPSILON_SIGMA_FICT = {
    "H": LennardJonesParams(2.2094300, 4.4778900, 0.5523570),
    "He": LennardJonesParams(1.9956100, 0.0009421, 0.4989030),
    "Li": LennardJonesParams(9.1228000, 1.0496900, 2.2807000),
    "Be": LennardJonesParams(6.8421000, 0.5729420, 1.7105300),
    "B": LennardJonesParams(6.0581100, 2.9670300, 1.5145300),
    "C": LennardJonesParams(5.4166600, 6.3695300, 1.3541700),
    "N": LennardJonesParams(5.0603000, 9.7537900, 1.2650800),
    "O": LennardJonesParams(4.7039500, 5.1264700, 1.1759900),
    "F": LennardJonesParams(4.0625000, 1.6059200, 1.0156200),
    "Ne": LennardJonesParams(4.1337700, 0.0036471, 1.0334400),
    "Na": LennardJonesParams(11.8311000, 0.7367450, 2.9577800),
    "Mg": LennardJonesParams(10.0493000, 0.0785788, 2.5123300),
    "Al": LennardJonesParams(8.6239000, 2.7006700, 2.1559700),
    "Si": LennardJonesParams(7.9111800, 3.1743100, 1.9778000),
    "P": LennardJonesParams(7.6260900, 5.0305000, 1.9065200),
    "S": LennardJonesParams(7.4835500, 4.3692700, 1.8708900),
    "Cl": LennardJonesParams(7.2697300, 4.4832800, 1.8174300),
    "Ar": LennardJonesParams(7.5548200, 0.0123529, 1.8887100),
    "K": LennardJonesParams(14.4682000, 0.5517990, 3.6170500),
    "Ca": LennardJonesParams(12.5439000, 0.1326790, 3.1359600),
    "Sc": LennardJonesParams(12.1162000, 1.6508000, 3.0290600),
    "Ti": LennardJonesParams(11.4035000, 1.1802700, 2.8508800),
    "V": LennardJonesParams(10.9046000, 2.7524900, 2.7261500),
    "Cr": LennardJonesParams(9.9067900, 1.5367900, 2.4767000),
    "Mn": LennardJonesParams(9.9067900, 0.5998880, 2.4767000),
    "Fe": LennardJonesParams(9.4078900, 1.1844200, 2.3519700),
    "Co": LennardJonesParams(8.9802600, 1.2776900, 2.2450600),
    "Ni": LennardJonesParams(8.8377200, 2.0757200, 2.2094300),
    "Cu": LennardJonesParams(9.4078900, 2.0446300, 2.3519700),
    "Zn": LennardJonesParams(8.6951700, 0.1915460, 2.1737900),
    "Ga": LennardJonesParams(8.6951700, 1.0642000, 2.1737900),
    "Ge": LennardJonesParams(8.5526300, 2.7017100, 2.1381600),
    "As": LennardJonesParams(8.4813600, 3.9599000, 2.1203400),
    "Se": LennardJonesParams(8.5526300, 3.3867700, 2.1381600),
    "Br": LennardJonesParams(8.5526300, 1.9706300, 2.1381600),
    "Kr": LennardJonesParams(8.2675400, 0.0173276, 2.0668900),
    "Rb": LennardJonesParams(15.6798000, 0.4682650, 3.9199500),
    "Sr": LennardJonesParams(13.8980000, 0.1339230, 3.4745100),
    "Y": LennardJonesParams(13.5417000, 2.7597500, 3.3854200),
    "Zr": LennardJonesParams(12.4726000, 3.0520100, 3.1181500),
    "Nb": LennardJonesParams(11.6886000, 5.2782000, 2.9221500),
    "Mo": LennardJonesParams(10.9759000, 4.4749900, 2.7439700),
    "Tc": LennardJonesParams(10.4770000, 3.3815900, 2.6192400),
    "Ru": LennardJonesParams(10.4057000, 1.9617200, 2.6014200),
    "Rh": LennardJonesParams(10.1206000, 2.4058200, 2.5301500),
    "Pd": LennardJonesParams(9.9067900, 1.3709700, 2.4767000),
    "Ag": LennardJonesParams(10.3344000, 1.6497600, 2.5836100),
    "Cd": LennardJonesParams(10.2632000, 0.0377447, 2.5657900),
    "In": LennardJonesParams(10.1206000, 0.8113140, 2.5301500),
    "Sn": LennardJonesParams(9.9067900, 1.9005700, 2.4767000),
    "Sb": LennardJonesParams(9.9067900, 3.0882800, 2.4767000),
    "Te": LennardJonesParams(9.8355200, 2.6312300, 2.4588800),
    "I": LennardJonesParams(9.9067900, 1.5393800, 2.4767000),
    "Xe": LennardJonesParams(9.9780700, 0.0238880, 2.4945200),
    "Cs": LennardJonesParams(17.3903000, 0.4166420, 4.3475900),
    "Ba": LennardJonesParams(15.3235000, 1.9000000, 3.8308600),
    "La": LennardJonesParams(14.7533000, 2.4996100, 3.6883200),
    "Ce": LennardJonesParams(14.5395000, 2.5700800, 3.6348700),
    "Pr": LennardJonesParams(14.4682000, 1.2994600, 3.6170500),
    "Nd": LennardJonesParams(14.3257000, 0.8196050, 3.5814100),
    "Pm": LennardJonesParams(14.1831000, 3.2413400, 3.5457800),
    "Sm": LennardJonesParams(14.1118000, 0.5211220, 3.5279600),
    "Eu": LennardJonesParams(14.1118000, 0.4299180, 3.5279600),
    "Gd": LennardJonesParams(13.9693000, 2.0995600, 3.4923200),
    "Tb": LennardJonesParams(13.8267000, 1.3999900, 3.4566900),
    "Dy": LennardJonesParams(13.6842000, 0.6900550, 3.4210500),
    "Ho": LennardJonesParams(13.6842000, 0.6900550, 3.4210500),
    "Er": LennardJonesParams(13.4704000, 0.7387660, 3.3676000),
    "Tm": LennardJonesParams(13.5417000, 0.5211220, 3.3854200),
    "Yb": LennardJonesParams(13.3278000, 0.1303990, 3.3319600),
    "Lu": LennardJonesParams(13.3278000, 1.4331500, 3.3319600),
    "Hf": LennardJonesParams(12.4726000, 3.3608600, 3.1181500),
    "Ta": LennardJonesParams(12.1162000, 4.0034300, 3.0290600),
    "W": LennardJonesParams(11.5460000, 6.8638900, 2.8865100),
    "Re": LennardJonesParams(10.7621000, 4.4387100, 2.6905100),
    "Os": LennardJonesParams(10.2632000, 4.2625300, 2.5657900),
    "Ir": LennardJonesParams(10.0493000, 3.7028700, 2.5123300),
    "Pt": LennardJonesParams(9.6929800, 3.1401000, 2.4232400),
    "Au": LennardJonesParams(9.6929800, 2.3058000, 2.4232400),
    "Hg": LennardJonesParams(9.4078900, 0.0454140, 2.3519700),
    "Tl": LennardJonesParams(10.3344000, 0.5770870, 2.5836100),
    "Pb": LennardJonesParams(10.4057000, 0.8589880, 2.6014200),
    "Bi": LennardJonesParams(10.5482000, 2.0798700, 2.6370600),
    "Po": LennardJonesParams(9.9780700, 1.8995300, 2.4945200),
    "At": LennardJonesParams(10.6908000, 1.3854420, 2.6727000),
    "Rn": LennardJonesParams(10.6908000, 0.0214992, 2.6727000),
    "Fr": LennardJonesParams(18.5307000, 0.3749778, 4.6326700),
    "Ra": LennardJonesParams(15.7511000, 1.7100000, 3.9377700),
    "Ac": LennardJonesParams(15.3235000, 2.2496490, 3.8308600),
    "Th": LennardJonesParams(14.6820000, 2.3130720, 3.6705000),
    "Pa": LennardJonesParams(14.2544000, 1.1695140, 3.5635900),
    "U": LennardJonesParams(13.9693000, 0.7376445, 3.4923200),
    "Np": LennardJonesParams(13.5417000, 2.9172060, 3.3854200),
    "Pu": LennardJonesParams(13.3278000, 0.4690098, 3.3319600),
    "Am": LennardJonesParams(12.8289000, 0.3869262, 3.2072400),
    "Cm": LennardJonesParams(12.0450000, 1.8896040, 3.0112400),
    "Bk": LennardJonesParams(11.9737000, 1.2599910, 2.9934200),
    "Cf": LennardJonesParams(11.9737000, 0.6210495, 2.9934200),
    "Es": LennardJonesParams(11.7599000, 0.6210495, 2.9399700),
    "Fm": LennardJonesParams(11.9024000, 0.6648894, 2.9756000),
    "Md": LennardJonesParams(12.3300000, 0.4690098, 3.0825100),
    "No": LennardJonesParams(12.5439000, 0.1173591, 3.1359600),
    "Lr": LennardJonesParams(11.4748000, 1.2898350, 2.8686900),
    "Rf": LennardJonesParams(11.1897000, 3.0247740, 2.7974200),
    "Db": LennardJonesParams(10.6195000, 3.6030870, 2.6548800),
    "Sg": LennardJonesParams(10.1919000, 6.1775010, 2.5479700),
    "Bh": LennardJonesParams(10.0493000, 3.9948390, 2.5123300),
    "Hs": LennardJonesParams(9.5504300, 3.8362770, 2.3876100),
    "Mt": LennardJonesParams(9.1940700, 3.3325830, 2.2985200),
    "Ds": LennardJonesParams(9.1228000, 2.8260900, 2.2807000),
    "Rg": LennardJonesParams(8.6239000, 2.0752200, 2.1559700),
    "Cn": LennardJonesParams(8.6951700, 0.0408726, 2.1737900),
    "Nh": LennardJonesParams(9.6929800, 0.5193783, 2.4232400),
    "Fl": LennardJonesParams(10.1919000, 0.7730892, 2.5479700),
    "Mc": LennardJonesParams(11.5460000, 1.8718830, 2.8865100),
    "Lv": LennardJonesParams(12.4726000, 1.7095770, 3.1181500),
    "Ts": LennardJonesParams(11.7599000, 1.2468978, 2.9399700),
    "Og": LennardJonesParams(11.1897000, 0.0193493, 2.7974200),
}


D_ELECTRONS_SPIN = [
    [0],
    [1],
    [2, 0],
    [3, 1],
    [4, 2, 0],
    [5, 3, 1],
    [4, 2, 0],
    [3, 1],
    [2, 0],
    [1],
    [0],
]


@dataclass
class f_d_s_p_electrons:
    f: int
    d: int
    s: int
    p: int


METAL_F_D_S_P_ELECTRONS = {
    "Li": f_d_s_p_electrons(0, 0, 1, 0),
    "Be": f_d_s_p_electrons(0, 0, 2, 0),
    "Na": f_d_s_p_electrons(0, 0, 1, 0),
    "Mg": f_d_s_p_electrons(0, 0, 2, 0),
    "Al": f_d_s_p_electrons(0, 0, 2, 1),
    "K": f_d_s_p_electrons(0, 0, 1, 0),
    "Ca": f_d_s_p_electrons(0, 0, 2, 0),
    "Sc": f_d_s_p_electrons(0, 1, 2, 0),
    "Ti": f_d_s_p_electrons(0, 2, 2, 0),
    "V": f_d_s_p_electrons(0, 3, 2, 0),
    "Cr": f_d_s_p_electrons(0, 5, 1, 0),
    "Mn": f_d_s_p_electrons(0, 5, 2, 0),
    "Fe": f_d_s_p_electrons(0, 6, 2, 0),
    "Co": f_d_s_p_electrons(0, 7, 2, 0),
    "Ni": f_d_s_p_electrons(0, 8, 2, 0),
    "Cu": f_d_s_p_electrons(0, 10, 1, 0),
    "Zn": f_d_s_p_electrons(0, 10, 2, 0),
    "Ga": f_d_s_p_electrons(0, 10, 2, 1),
    "Ge": f_d_s_p_electrons(0, 10, 2, 2),
    "Rb": f_d_s_p_electrons(0, 0, 1, 0),
    "Sr": f_d_s_p_electrons(0, 0, 2, 0),
    "Y": f_d_s_p_electrons(0, 1, 2, 0),
    "Zr": f_d_s_p_electrons(0, 2, 2, 0),
    "Nb": f_d_s_p_electrons(0, 4, 1, 0),
    "Mo": f_d_s_p_electrons(0, 5, 1, 0),
    "Tc": f_d_s_p_electrons(0, 5, 2, 0),
    "Ru": f_d_s_p_electrons(0, 7, 1, 0),
    "Rh": f_d_s_p_electrons(0, 8, 1, 0),
    "Pd": f_d_s_p_electrons(0, 10, 0, 0),
    "Ag": f_d_s_p_electrons(0, 10, 1, 0),
    "Cd": f_d_s_p_electrons(0, 10, 2, 0),
    "In": f_d_s_p_electrons(0, 10, 2, 1),
    "Sn": f_d_s_p_electrons(0, 10, 2, 2),
    "Sb": f_d_s_p_electrons(0, 10, 2, 3),
    "Cs": f_d_s_p_electrons(0, 0, 1, 0),
    "Ba": f_d_s_p_electrons(0, 0, 2, 0),
    "La": f_d_s_p_electrons(0, 1, 2, 0),
    "Ce": f_d_s_p_electrons(1, 1, 2, 0),
    "Pr": f_d_s_p_electrons(3, 0, 2, 0),
    "Nd": f_d_s_p_electrons(4, 0, 2, 0),
    "Pm": f_d_s_p_electrons(5, 0, 2, 0),
    "Sm": f_d_s_p_electrons(6, 0, 2, 0),
    "Eu": f_d_s_p_electrons(7, 0, 2, 0),
    "Gd": f_d_s_p_electrons(7, 1, 2, 0),
    "Tb": f_d_s_p_electrons(9, 0, 2, 0),
    "Dy": f_d_s_p_electrons(10, 0, 2, 0),
    "Ho": f_d_s_p_electrons(11, 0, 2, 0),
    "Er": f_d_s_p_electrons(12, 0, 2, 0),
    "Tm": f_d_s_p_electrons(13, 0, 2, 0),
    "Yb": f_d_s_p_electrons(14, 0, 2, 0),
    "Lu": f_d_s_p_electrons(14, 1, 2, 0),
    "Hf": f_d_s_p_electrons(14, 2, 2, 0),
    "Ta": f_d_s_p_electrons(14, 3, 2, 0),
    "W": f_d_s_p_electrons(14, 4, 2, 0),
    "Re": f_d_s_p_electrons(14, 5, 2, 0),
    "Os": f_d_s_p_electrons(14, 6, 2, 0),
    "Ir": f_d_s_p_electrons(14, 7, 2, 0),
    "Pt": f_d_s_p_electrons(14, 9, 1, 0),
    "Au": f_d_s_p_electrons(14, 10, 1, 0),
    "Hg": f_d_s_p_electrons(14, 10, 2, 0),
    "Tl": f_d_s_p_electrons(14, 10, 2, 1),
    "Pb": f_d_s_p_electrons(14, 10, 2, 2),
    "Bi": f_d_s_p_electrons(14, 10, 2, 3),
}

METAL_F_D_S_P_ELECTRONS_KEYS = list(METAL_F_D_S_P_ELECTRONS.keys())

METAL_VALENCE_AVAILABLE_PRIOR = {
    "Li": [1],
    "Na": [1],
    "K": [1],
    "Rb": [1],
    "Cs": [1],
    "Be": [2],
    "Mg": [2],
    "Ca": [2],
    "Sr": [2],
    "Ba": [2],
    "Sc": [3],
    "Y": [3],
    "La": [3],
    "Ce": [3],
    "Pr": [3],
    "Nd": [3],
    "Pm": [3],
    "Sm": [3],
    "Eu": [3],
    "Gd": [3],
    "Tb": [3],
    "Dy": [3],
    "Ho": [3],
    "Er": [3],
    "Tm": [3],
    "Yb": [3],
    "Lu": [3],
    "Ti": [4, 3],
    "Zr": [4],
    "Hf": [4],
    "V": [5, 4, 3, 2],
    "Nb": [5],
    "Ta": [5],
    "Cr": [6, 3, 2],
    "Mo": [3, 4, 6],
    "W": [6, 4],
    "Mn": [7, 6, 4, 3, 2],
    "Tc": [7, 4],
    "Re": [7, 4],
    "Fe": [2, 3, 4, 0],
    "Ru": [4, 3, 2, 0],
    "Os": [8, 6, 4, 0],
    "Co": [3, 2, 0],
    "Rh": [3, 0],
    "Ir": [4, 3, 0],
    "Ni": [2, 0],
    "Pd": [2, 4, 0],
    "Pt": [4, 2, 0],
    "Cu": [2, 1, 3],
    "Ag": [1],
    "Au": [3, 1],
    "Zn": [2],
    "Cd": [2],
    "Hg": [2, 1],
    "Al": [3],
    "Ga": [3],
    "In": [3],
    "Tl": [3],
    "Ge": [4],
    "Sn": [4],
    "Pb": [4],
    "Sb": [5, 3],
    "Bi": [5, 3],
}

METAL_VALENCE_AVAILABLE_MINOR = {
    "Li": [0],
    "Na": [0],
    "K": [0],
    "Rb": [0],
    "Cs": [0],
    "Be": [0],
    "Mg": [0],
    "Ca": [0],
    "Sr": [0],
    "Ba": [0],
    "Sc": [2, 1, 0],
    "Y": [2, 1, 0],
    "La": [2, 1, 0],
    "Ce": [2, 1, 0],
    "Pr": [2, 1, 0],
    "Nd": [2, 1, 0],
    "Pm": [2, 1, 0],
    "Sm": [2, 1, 0],
    "Eu": [2, 1, 0],
    "Gd": [2, 1, 0],
    "Tb": [2, 1, 0],
    "Dy": [2, 1, 0],
    "Ho": [2, 1, 0],
    "Er": [2, 1, 0],
    "Tm": [2, 1, 0],
    "Yb": [2, 1, 0],
    "Lu": [2, 1, 0],
    "Ti": [1, -1, -2, 0],
    "Zr": [3, 2, 1, -2, 0],
    "Hf": [3, 2, 1, -2, 0],
    "V": [1, -1, -3, 0],
    "Nb": [4, 3, 2, 1, -1, -3, 0],
    "Ta": [4, 3, 2, 1, -1, -3, 0],
    "Cr": [5, 4, 1, -1, -2, -4, 0],
    "Mo": [5, 2, 1, -1, -2, -4, 0],
    "W": [5, 3, 2, 1, -1, -2, -4, 0],
    "Mn": [0, 5, 1, -1, -2, -3, 0],
    "Tc": [6, 3, 2, 5, 1, -1, -3, 0],
    "Re": [6, 3, 2, 5, 1, -1, -3, 0],
    "Fe": [6, 7, 5, 1, -1, -2, -4],
    "Ru": [2, 8, 7, 6, 5, 1, -2, -4],
    "Os": [3, 2, 7, 5, 1, -1, -2, -4],
    "Co": [5, 4, 1, -1, -3],
    "Rh": [2, 7, 6, 5, 4, 1, -1, -3],
    "Ir": [2, 9, 8, 7, 6, 5, 1, -1, -3],
    "Ni": [3, 4, 1, -1, -2],
    "Pd": [0, 5, 3, 1],
    "Pt": [0, 6, 5, 3, 1, -1, -2, 0],
    "Cu": [4, 0],
    "Ag": [2, 4, 3, 0],
    "Au": [5, 2, -1, 0],
    "Zn": [1, 0],
    "Cd": [1, 0],
    "Hg": [0],
    "Al": [0],
    "Ga": [0],
    "In": [0],
    "Tl": [0],
    "Ge": [0],
    "Sn": [0],
    "Pb": [0],
    "Sb": [0],
    "Bi": [0],
}


def get_possible_metal_radicals(metal: str, valence: int) -> set[int]:
    """
    Get possible radicals for a given metal and valence.

    If the valence is less than or equal to the sum of the s and p electrons,
    then the possible radicals are the set of electrons with the same spin as the
    valence minus the sum of the s and p electrons.

    If the valence is less than or equal to the sum of the s, p, and d electrons,
    then the possible radicals are the set of electrons with the same spin as the
    valence minus the sum of the s, p, and d electrons.

    If the valence is less than or equal to the sum of the s, p, d, and f electrons,
    then the possible radicals are the set of electrons with the same spin as the
    valence minus the sum of the s, p, d, and f electrons.

    Parameters:

        metal (str): The metal symbol.
        valence (int): The valence of the metal.

    Returns:
        set[int]: A set of possible radicals.
    """
    f_d_s_p = METAL_F_D_S_P_ELECTRONS[metal]
    f, d, s, p = f_d_s_p.f, f_d_s_p.d, f_d_s_p.s, f_d_s_p.p
    if valence <= s + p:
        return {(f + s + p - valence) % 2 + dd for dd in D_ELECTRONS_SPIN[d]}
    if valence <= s + p + d:
        return {f % 2 + dd for dd in D_ELECTRONS_SPIN[d - valence + s + p]}
    if valence <= s + p + d + f:
        return {f % 2}
    raise ValueError("Valence is too high for this metal")


@dataclass
class element:
    symbol: str
    name: str
    atomic_number: int
    atomic_mass: float
    num_outer_electrons: int
    default_valence: int


NON_METAL_DICT = {
    1: element("H", "Hydrogen", 1, 1.008, 1, 1),
    2: element("He", "Helium", 2, 4.003, 2, 0),
    5: element("B", "Boron", 5, 10.812, 3, 3),
    6: element("C", "Carbon", 6, 12.011, 4, 4),
    7: element("N", "Nitrogen", 7, 14.007, 5, 3),
    8: element("O", "Oxygen", 8, 15.999, 6, 2),
    9: element("F", "Fluorine", 9, 18.998, 7, 1),
    10: element("Ne", "Neon", 10, 20.18, 8, 0),
    14: element("Si", "Silicon", 14, 28.086, 4, 4),
    15: element("P", "Phosphorus", 15, 30.974, 5, 3),
    16: element("S", "Sulfur", 16, 32.066, 6, 2),
    17: element("Cl", "Chlorine", 17, 35.453, 7, 1),
    18: element("Ar", "Argon", 18, 39.948, 8, 0),
    32: element("Ge", "Germanium", 32, 72.61, 4, 4),
    33: element("As", "Arsenic", 33, 74.922, 5, 3),
    34: element("Se", "Selenium", 34, 78.96, 6, 2),
    35: element("Br", "Bromine", 35, 79.904, 7, 1),
    36: element("Kr", "Krypton", 36, 83.798, 8, 0),
    51: element("Sb", "Antimony", 51, 121.76, 5, 3),
    52: element("Te", "Tellurium", 52, 127.6, 6, 2),
    53: element("I", "Iodine", 53, 126.904, 7, 1),
    54: element("Xe", "Xenon", 54, 131.29, 8, 0),
    84: element("Po", "Polonium", 84, 209.0, 6, 2),
    85: element("At", "Astatine", 85, 210.0, 7, 1),
    86: element("Rn", "Radon", 86, 222.0, 8, 0),
}

HETEROATOM = (9, 8, 17, 7, 35, 54, 16, 34, 15)
