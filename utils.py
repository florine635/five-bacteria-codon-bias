# -*- coding: utf-8 -*-
"""
utils.py —— 公共函数：密码子表、CDS 清洗、三个核心指标
=====================================================

三个指标的定义（报告"方法"一节直接照抄这三段即可）：

1) RSCU（Relative Synonymous Codon Usage，相对同义密码子使用度）
        RSCU(密码子 i) = X_i / ( (1/k) * Σ_j X_j )
   其中 k 是同一个氨基酸对应的同义密码子个数，X_i 是该密码子的观测次数。
   - RSCU = 1  → 该密码子被随机使用
   - RSCU > 1  → 偏好使用
   - RSCU < 1  → 回避使用
   只对简并度为 2/3/4/6 的氨基酸有意义；Met、Trp 只有唯一密码子，不计算。

2) GC3s（第三位密码子 GC 含量）
   只看第 3 位碱基，统计 G/C 的比例。排除 Met、Trp 和终止密码子。
   它衡量"突变压力"这一侧的影响。

3) ENC（Effective Number of Codons，有效密码子数，Wright 1990）
        F_k = ( n * Σ p_i² - 1 ) / ( n - 1 )       # 某简并度 k 家族内的纯合度
        ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
   其中 n 是该家族的密码子总数，p_i 是各同义密码子的频率，
   F2/F3/F4/F6 是在对应简并度家族内取平均。
   ENC 取值范围 20（极端偏好，每个氨基酸只用 1 个密码子）到 61（完全随机）。
   系数 9/1/5/3 正好是标准遗传密码表里 2/3/4/6 简并家族的个数。

4) ENC 期望曲线（画 ENC-plot 用，Wright 1990）
        ENC_exp = 2 + GC3s + 29 / ( GC3s² + (1 - GC3s)² )
   这是一条"仅由突变压力决定、没有选择"的曲线。实际点落在曲线下方，
   说明该基因受到了翻译选择（偏好高表达所需密码子）。
"""

from collections import Counter, defaultdict

# ===========================================================================
# 遗传密码表：优先用 Biopython；没装 Biopython 时用内置表
# ===========================================================================
# 为什么要有内置表：
#   04_plot.py 和 05_gui.py 只需要"每个密码子属于哪个氨基酸"这一张对照表，
#   既不联网、也不解析序列。如果为了这张表强制安装 Biopython，
#   会让"只想打开界面看图"的人也被卡住（这是本项目最初的真实 bug）。
#   所以只有 01_fetch_data.py（联网取数）才真正依赖 Biopython。
#
# 内置表怎么保证没写错：
#   用 NCBI 遗传密码表的经典四行串（三行碱基 + 一行氨基酸）构造，
#   并在 verify_codon_table() 里与 Biopython 逐密码子对拍。
#   装了 Biopython 的机器上跑一次 utils.py 就会做这个校验。

try:
    from Bio.Data import CodonTable as _bio_codontable
    HAS_BIOPYTHON = True
except Exception:                                     # noqa: BLE001
    _bio_codontable = None
    HAS_BIOPYTHON = False

# NCBI 表 11（细菌、古菌与质体）的规范写法，顺序不能变
_NCBI_BASE1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_NCBI_BASE2 = "TTTTCCCCAAAAGGGG" * 4
_NCBI_BASE3 = "TCAG" * 16
_NCBI_AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

# 允许的碱基
_ACGT = set("ACGT")


class _MiniCodonTable:
    """最小密码表对象，接口与 Bio.Data.CodonTable 需要的那两个属性一致。"""

    def __init__(self, forward_table, stop_codons):
        self.forward_table = forward_table      # 密码子 -> 氨基酸单字母
        self.stop_codons = stop_codons          # 终止密码子列表

    def __repr__(self):
        return (f"<内置密码表：{len(self.forward_table)} 个有义密码子 + "
                f"{len(self.stop_codons)} 个终止密码子>")


def _build_builtin_table():
    fwd, stops = {}, []
    for b1, b2, b3, aa in zip(_NCBI_BASE1, _NCBI_BASE2, _NCBI_BASE3, _NCBI_AA_11):
        codon = b1 + b2 + b3
        if aa == "*":
            stops.append(codon)
        else:
            fwd[codon] = aa
    return _MiniCodonTable(fwd, stops)


_BUILTIN_TABLE = _build_builtin_table()


def get_table(table_id=11):
    """
    按 NCBI 编号取遗传密码表。

    装了 Biopython 就返回 Biopython 的表（支持全部编号）；
    没装则返回内置的标准表，此时 table_id 只支持 1 和 11
    （这两个编号在有义/终止密码子上完全一致，本项目 5 个物种都用它）。
    """
    if HAS_BIOPYTHON:
        return _bio_codontable.unambiguous_dna_by_id[table_id]
    if table_id not in (1, 11):
        raise ValueError(
            f"没装 Biopython 时只能用内置密码表（table_id = 1 或 11），"
            f"当前传入 {table_id}。\n"
            f"需要其它遗传密码表（如支原体用 4、线粒体用 2）请先安装：\n"
            f"    pip install -i https://pypi.tuna.tsinghua.edu.cn/simple biopython")
    return _BUILTIN_TABLE


def verify_codon_table():
    """
    把内置表与 Biopython 的表逐密码子对拍。

    返回 (是否一致, 说明文字)。没装 Biopython 时返回 (None, 提示)。
    这是防止"手写密码表打错一个字"的唯一可靠手段 —— 一旦打错，
    RSCU/ENC 会整体算错，而且不会报错。
    """
    if not HAS_BIOPYTHON:
        return None, "未安装 Biopython，无法对拍（内置表仍在工作）"

    ref = _bio_codontable.unambiguous_dna_by_id[11]
    diff_fwd = {c: (ref.forward_table.get(c), _BUILTIN_TABLE.forward_table.get(c))
                for c in set(ref.forward_table) | set(_BUILTIN_TABLE.forward_table)
                if ref.forward_table.get(c) != _BUILTIN_TABLE.forward_table.get(c)}
    diff_stop = set(ref.stop_codons) ^ set(_BUILTIN_TABLE.stop_codons)
    if diff_fwd or diff_stop:
        return False, f"不一致！有义密码子差异 {diff_fwd}，终止密码子差异 {diff_stop}"
    return True, (f"一致：{len(ref.forward_table)} 个有义密码子 + "
                  f"{len(ref.stop_codons)} 个终止密码子全部命中")


def codon_families(table):
    """
    返回 (families, stops)。

    families: dict 氨基酸 -> 该氨基酸的同义密码子列表
              （只保留简并度 >= 2 的家族，Met/Trp 这类单密码子氨基酸排除）
    stops   : 终止密码子集合
    """
    fam = defaultdict(list)
    for codon, aa in table.forward_table.items():
        fam[aa].append(codon)
    fam = {aa: sorted(v) for aa, v in fam.items() if len(v) >= 2}
    return fam, set(table.stop_codons)


def count_codons(codons):
    """密码子列表 -> Counter。"""
    return Counter(codons)


def gc3s(codons, table):
    """
    第三位密码子 GC 含量。codons 不含终止密码子。

    按惯例排除 Met(ATG) 与 Trp(TGG)：它们只有一个密码子，第 3 位不含信息。
    """
    keep = [c for c in codons
            if c in table.forward_table and table.forward_table[c] not in ("M", "W")]
    if not keep:
        return float("nan")
    return sum(1 for c in keep if c[2] in "GC") / len(keep)


def rscu_vector(counts, families):
    """
    计算 59 个同义密码子的 RSCU，返回 dict 密码子 -> 值。

    counts 是 Counter；families 来自 codon_families()。
    某家族计数为 0 时整族跳过（该基因没有这个氨基酸）。
    """
    out = {}
    for aa, codons in families.items():
        total = sum(counts.get(c, 0) for c in codons)
        if total == 0:
            continue
        k = len(codons)
        expected = total / k
        for c in codons:
            out[c] = counts.get(c, 0) / expected
    return out


# ENC 的理论取值区间（20 = 极端偏好，61 = 完全随机）
ENC_MIN, ENC_MAX = 20.0, 61.0


def enc_value(counts, families):
    """
    ENC：有效密码子数（Wright 1990）。

    家族内纯合度 F_k 越小 → 密码子使用越均匀 → 该家族贡献的 ENC 越大。

    ⚠️ 返回值**不做截断**，可能超过理论上限 61。这不是 bug，是这个估计量的已知缺陷：
       当某个高简并家族（6 简并的 Leu/Arg/Ser）在该基因里只有极少数密码子时，
       F_k 会趋近 0，而 ENC 公式里要除以 F_k（3/F6、9/F2），结果会被放大到几百。
       正确做法是**把这类基因剔除**，而不是把数字按到 61 —— 02_compute_metrics.py
       会做这个判断，并记进过滤日志（原因 enc_out_of_range）。
       报告里要写清楚剔除了多少条、为什么，否则别人复现不出你的数字。
    """
    fk = defaultdict(list)          # 简并度 -> [F_k, ...]
    for aa, codons in families.items():
        total = sum(counts.get(c, 0) for c in codons)
        k = len(codons)
        if total <= 1:              # n<=1 时 F 无定义，跳过
            continue
        p2 = sum((counts.get(c, 0) / total) ** 2 for c in codons)
        fk[k].append((total * p2 - 1) / (total - 1))

    # 每个简并度取平均，然后按 Wright 公式累加
    mean_f = {}
    for k, vals in fk.items():
        if vals:
            mean_f[k] = sum(vals) / len(vals)

    enc = 2.0
    for k, coef in ((2, 9), (3, 1), (4, 5), (6, 3)):
        f = mean_f.get(k)
        if f and f > 0:
            enc += coef / f
    return enc


def enc_expected(gc3s_array):
    """ENC 期望曲线（仅突变压力，无选择）。输入可以是 float 或 numpy 数组。"""
    return 2 + gc3s_array + 29.0 / (gc3s_array ** 2 + (1 - gc3s_array) ** 2)


def clean_cds(seq, table, min_codons=100, require_atg=True, min_acgt=0.99):
    """
    清洗一条 CDS，返回 (codons, reason)。

    codons : 去掉终止密码子后的密码子列表；不合格时返回 None
    reason : 'ok' 或不合格原因，用于过滤日志

    过滤规则（报告"方法"里要写清楚，否则结果不可复现）：
        - 模糊碱基过多 → 'ambiguous'
        - 长度不是 3 的倍数 → 'not_multiple_of_3'
        - 起始不是 ATG（可关闭）→ 'not_atg_start'
        - 末尾不是终止密码子 → 'no_stop'
        - 中间有终止密码子 → 'internal_stop'
        - 密码子数太少 → 'too_short'
    """
    seq = seq.upper().replace("U", "T")
    if not seq:
        return None, "empty"
    if sum(1 for b in seq if b in _ACGT) / len(seq) < min_acgt:
        return None, "ambiguous"
    if len(seq) % 3 != 0:
        return None, "not_multiple_of_3"
    if require_atg and not seq.startswith("ATG"):
        return None, "not_atg_start"

    stops = set(table.stop_codons)
    if seq[-3:] not in stops:
        return None, "no_stop"

    body = [seq[i:i + 3] for i in range(0, len(seq) - 3, 3)]
    if any(c in stops for c in body):
        return None, "internal_stop"
    if len(body) < min_codons:
        return None, "too_short"
    return body, "ok"


def read_fasta(path):
    """
    极简 FASTA 解析（不依赖 Biopython，方便在其他脚本里复用）。

    返回 (header_without_gt, sequence) 的生成器。
    """
    header, chunks = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def parse_gene(header):
    """从 fasta_cds_na 的 header 里抠出 [gene=xxx]，取不到返回空串。"""
    for tok in header.split():
        if tok.startswith("[gene="):
            return tok[6:].rstrip("]")
    return ""


def setup_chinese_font():
    """
    让 matplotlib 能显示中文（否则标签全是方框）。
    返回实际生效的字体名，找不到就返回 None。

    Windows 自带 SimHei / 微软雅黑；macOS 用 PingFang；Linux 通常要自己装。
    """
    import matplotlib
    from matplotlib import font_manager

    candidates = ["SimHei", "Microsoft YaHei", "PingFang SC",
                  "Noto Sans CJK SC", "WenQuanYi Zen Hei", "Arial Unicode MS"]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            matplotlib.rcParams["font.sans-serif"] = [name] + \
                matplotlib.rcParams["font.sans-serif"]
            matplotlib.rcParams["axes.unicode_minus"] = False   # 负号别变方框
            return name
    return None


if __name__ == "__main__":
    # 自检 1：内置密码表与 Biopython 是否逐项一致（防手打错字）
    ok, msg = verify_codon_table()
    print(f"[密码表对拍] {msg}")
    print(f"[环境] Biopython {'已安装' if HAS_BIOPYTHON else '未安装（用内置表）'}")
    if ok is False:
        raise SystemExit("内置密码表与 Biopython 不一致，禁止继续！")

    # 自检 2：用一段已知序列验证三个指标的实现没跑偏
    t = get_table(11)
    fam, stops = codon_families(t)
    print("简并家族数:", len(fam), " 终止密码子:", sorted(stops))

    demo = ("ATG" + "GCT" * 20 + "GAA" * 20 + "CTG" * 20 + "TGG" + "TAA")
    body, reason = clean_cds(demo, t, min_codons=10)
    c = count_codons(body)
    print("clean_cds ->", reason, " 密码子数:", len(body))
    print("ENC =", round(enc_value(c, fam), 3), " GC3s =", round(gc3s(body, t), 3))
    print("理论极值：全同义密码子各用一次时 ENC≈61，极端偏好时 ENC≈20")
