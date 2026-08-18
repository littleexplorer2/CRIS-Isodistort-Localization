"""
中英双语科学术语对照表（Terminology Dictionary）

英文术语对齐优先级：
1. ISODISTORT 官网（https://iso.byu.edu/isodistorthelp.php）用语
2. 《International Tables for Crystallography》（ITC）英文原版

中文术语对齐优先级：
1. 《晶体学名词》（全国科学技术名词审定委员会审定，1996）
2. 《国际晶体学表》（ITC）中文译本

每个词条附带来源标注：
- [web]    isodistort 官网帮助页
- [ITC]    International Tables for Crystallography
- [CT]     《晶体学名词》
- [com]    领域通用（相变/群论文献通行用法，无单一权威出处）

使用方式：
    from isocore.i18n import term_en2zh, term_zh2en

注意：本表只收录“科学术语”本身，不含界面提示文案（界面文案见 messages.py）。
"""

# 英文 → 中文（以《晶体学名词》为准；未收录者用通行译法并注明）
TERMS_EN2ZH: dict[str, str] = {
    # ---- 晶体学基础（CT/ITC）----
    "space group": "空间群",                       # CT
    "point group": "点群",                         # CT
    "crystal system": "晶系",                      # CT
    "triclinic": "三斜",                           # CT
    "monoclinic": "单斜",                          # CT
    "orthorhombic": "正交",                        # CT
    "tetragonal": "四方",                          # CT
    "trigonal": "三方",                            # CT
    "hexagonal": "六方",                           # CT
    "cubic": "立方",                               # CT
    "lattice": "点阵",                             # CT（“晶格”为通行同义词）
    "unit cell": "晶胞",                           # CT
    "primitive cell": "原胞",                      # CT
    "conventional cell": "惯用晶胞",               # CT
    "supercell": "超胞",                           # [web] 官网超胞输入框用语
    "superlattice": "超晶格",                      # CT
    "sublattice": "子格",                          # [web]
    "superstructure": "超结构",                    # CT
    "fractional coordinates": "分数坐标",          # CT
    "lattice parameters": "晶格参数",              # CT
    "lattice constant": "晶格常数",                # CT
    "interatomic distance": "原子间距",            # [web]
    "nearest-neighbor distance": "最近邻距离",     # [web]
    "crystallographic axis": "晶轴",               # CT

    # ---- 空间群 / 对称性（ITC）----
    "space group number": "空间群序号",            # [web]/[ITC]
    "space group symbol": "空间群符号",            # ITC
    "Hermann-Mauguin symbol": "赫尔曼-莫甘符号",   # CT
    "site symmetry": "位置对称性",                 # CT
    "general position": "一般位置",                # CT
    "special position": "特殊位置",                # CT
    "Wyckoff position": "Wyckoff 位置",            # CT（附录表）
    "maximal subgroup": "极大子群",                # CT
    "minimal supergroup": "极小子群",              # CT
    "subgroup index": "子群指数",                  # [ITC]/[com]
    "origin choice": "原点选择",                   # ITC
    "cell choice": "胞型选择",                     # ITC
    "setting": "取位",                             # ITC（空间群“取位”）
    "unique axis": "特征轴",                       # CT（单斜唯一轴）
    "centering": "带心",                           # CT（如 C 带心/F 带心）
    "basis vectors": "基矢",                       # ITC
    "reciprocal lattice": "倒易点阵",              # CT
    "Brillouin zone": "布里渊区",                  # CT

    # ---- 相变 / 畸变（[web] + 通行物理术语）----
    "parent structure": "母体结构",                # [web]
    "parent phase": "母相",                        # [com] 相变文献
    "distorted structure": "畸变结构",             # [web]
    "distortion": "畸变",                          # [web]
    "displacement": "位移",                        # [web]
    "displacive": "位移型",                        # [web]（官网 Displacive 畸变）
    "atomic displacement": "原子位移",             # [web]
    "order": "有序化",                             # [web]（原子有序）
    "occupational": "占据型",                      # [web]（官网 Occupational 畸变）
    "rotational": "转动型",                        # [web]（官网 Rotational 畸变）
    "atomic ordering": "原子有序化",               # [web]
    "strain": "应变",                              # [web]
    "lattice strain": "晶格应变",                  # [web]
    "magnetic moment": "磁矩",                     # [web]/[CT]
    "all": "全部",                                 # [web]（作用域复选框）
    "none": "无",                                  # [web]（作用域复选框）
    "mode": "模式",                                # [web]
    "symmetry mode": "对称模式",                   # [web]
    "mode amplitude": "模式幅度",                  # [web]
    "amplitude": "幅度",                           # [com] 结构精修/相变文献
    "subgroup": "子群",                            # [ITC]
    "order parameter": "序参量",                   # [com] 群论物理标准
    "order parameter direction (OPD)": "序参量方向 (OPD)",  # [web]
    "isotropy subgroup": "各向同性子群",           # [web]/Stokes-Hatch
    "primary order parameter": "主序参量",         # [web]
    "secondary order parameter": "次级序参量",     # [web]
    "superposed order parameters": "叠加序参量",   # [web]
    "independent modulation": "独立调制",          # [web]
    "modulation": "调制",                          # [web]
    "commensurate": "公度",                        # [web]/[ITC]
    "incommensurate": "无公度",                    # [web]/[ITC]
    "superspace group": "超空间群",                # [web]/[ITC]
    "superspace": "超空间",                        # [web]/[ITC]
    "domain": "畴",                                # CT
    "domain generator": "畴生成元",                # [web]
    "domain count": "畴数",                        # [com]
    "origin shift": "原点平移",                    # [ITC]
    "lattice orientation": "点阵取向",             # [web]
    "internal orientation": "内部取向",            # [web]
    "ferroelectric": "铁电",                       # [web]
    "ferroelastic": "铁弹",                        # [web]
    "polar mode": "极性模式",                      # [com]
    "octahedral tilt": "八面体倾转",               # [com] 钙钛矿文献
    "mode decomposition": "模式分解",              # [web]
    "atom matching method": "原子匹配方法",        # [web]
    "nearest-site method": "最近位置法",           # [web]
    "robust method": "稳健法",                     # [web]
    "decomposition residual": "分解残差",          # [com]

    # ---- 不可约表示 / k 点（[web]/[ITC]/[CT]）----
    "irreducible representation (IR)": "不可约表示 (IR)",  # CT
    "irrep": "不可约表示",                         # [web]
    "k vector": "k 矢量",                          # [web]/[ITC]
    "k point": "k 点",                             # [web]
    "special k point": "特殊 k 点",                # [web]
    "star of k": "k 的星",                         # [ITC]
    "arm of the star": "星臂",                     # [ITC]
    "representation space": "表示空间",            # [com]
    "dimension": "维数",                           # CT
    "degeneracy": "简并度",                        # [ITC]
    "translational mode": "平移模式",              # [web]
    "longitudinal mode": "纵模式",                 # [ITC]
    "transverse mode": "横模式",                   # [ITC]
    "IR active": "红外活性",                       # [ITC]
    "Raman active": "拉曼活性",                    # [ITC]
    "Landau theory": "朗道理论",                   # [com]
    "renormalization group": "重正化群",           # [com]

    # ---- 物理畸变类型（[web] 复选框）----
    "displacement distortion": "位移畸变",         # [web]
    "occupancy distortion": "占据率畸变",          # [web]
    "strain distortion": "应变畸变",               # [web]
    "magnetic distortion": "磁性畸变",             # [web]
    "magnetic": "磁性",                            # CT
    "occupancy": "占据率",                         # CT

    # ---- 结构输出 / 文件（[web]）----
    "CIF file": "CIF 文件",                        # [web]
    "distortion file": "畸变文件",                 # [web]
    "interactive visualization": "交互可视化",     # [web]
    "subgroup tree": "子群树",                     # [web]（Bärnighausen 树）
    "refinement": "精修",                          # [com] 结构精修
}

# 中文 → 英文（自动由英文表反转生成）
TERMS_ZH2EN: dict[str, str] = {v: k for k, v in TERMS_EN2ZH.items()}


def term_en2zh(english: str) -> str:
    """英文科学术语 → 中文（《晶体学名词》优先）。未收录时原样返回。"""
    return TERMS_EN2ZH.get(english, english)


def term_zh2en(chinese: str) -> str:
    """中文科学术语 → 英文（用于切换回英文展示）。未收录时原样返回。"""
    return TERMS_ZH2EN.get(chinese, chinese)


def translate_term(term: str, language: str) -> str:
    """按当前语言翻译科学术语：en 返回英文原词，zh 返回中文。"""
    if language == "zh":
        return term_en2zh(term)
    return term
