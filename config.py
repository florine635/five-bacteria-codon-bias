# -*- coding: utf-8 -*-
"""
config.py —— 全局配置（改参数只改这里）
=====================================

换物种、改过滤阈值、改配色，全部集中在本文件。

⚠️ 关于登录号
    下面 5 个染色体登录号已经用 NCBI esummary 核实过（物种名、长度都对得上）。
    你自己新增物种时，**不要凭记忆填登录号**，先跑：

        python 01_fetch_data.py --verify

    它会把每个登录号的物种名和长度打印出来给你核对。
"""

from pathlib import Path

# ---------------------------------------------------------------- 目录
ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"          # 原始 FASTA 缓存（联网取一次，之后离线可跑）
FIG_DIR = ROOT / "figures"        # 输出图片
OUT_DIR = ROOT / "output"         # 输出表格 / 中间结果

for _d in (DATA_DIR, FIG_DIR, OUT_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- NCBI
NCBI_EMAIL = "your_name@example.edu"   # ← 改成你自己的邮箱，NCBI 要求提供
NCBI_RETRY = 3                          # 每个请求重试次数（限流/连接重置时很有用）
NCBI_SLEEP = 0.4                        # 请求间隔（秒），别调太小

# ---------------------------------------------------------------- 物种
# key = 图例名；acc = 染色体登录号；table = NCBI 遗传密码表编号
# gc_label 只是给人看的 GC 档位说明，不参与计算
SPECIES = {
    "B. subtilis":     {"acc": "NC_000964.3", "table": 11, "gc_label": "低 GC"},
    "E. coli":         {"acc": "NC_000913.3", "table": 11, "gc_label": "中 GC"},
    "M. tuberculosis": {"acc": "NC_000962.3", "table": 11, "gc_label": "高 GC"},
    "P. aeruginosa":   {"acc": "NC_002516.2", "table": 11, "gc_label": "高 GC"},
    "T. thermophilus": {"acc": "NC_006461.1", "table": 11, "gc_label": "极高 GC"},
}

# 展示顺序：按 GC3s 从低到高，图里更好读
SPECIES_ORDER = ["B. subtilis", "E. coli", "M. tuberculosis",
                 "P. aeruginosa", "T. thermophilus"]

# Okabe-Ito 色盲友好配色（不会红绿对撞）
COLORS = {
    "B. subtilis":     "#0072B2",
    "E. coli":         "#009E73",
    "M. tuberculosis": "#D55E00",
    "P. aeruginosa":   "#CC79A7",
    "T. thermophilus": "#E69F00",
}

# ---------------------------------------------------------------- CDS 过滤阈值
MIN_CODONS = 100        # 少于这么多密码子的 CDS 丢掉（太短，RSCU/ENC 不可靠）
REQUIRE_ATG = True      # 只保留 ATG 起始（见报告"局限性"：会漏掉 GTG/TTG 起始基因）
MIN_FRACTION_ACGT = 0.99  # 含过多 N 等模糊碱基的序列丢掉

# ---------------------------------------------------------------- 输出文件
CSV_CDS = OUT_DIR / "codon_metrics_per_CDS.csv"     # 每条 CDS 一行
CSV_SUMMARY = OUT_DIR / "species_summary.csv"       # 每个物种一行
CSV_FILTERLOG = OUT_DIR / "filter_log.csv"          # 过滤统计（报告要写）
CSV_RSCU = OUT_DIR / "rscu_matrix.csv"              # 物种 × 密码子 的 RSCU 均值
CSV_PCA = OUT_DIR / "pca_coordinates.csv"           # PCA 坐标
JSON_ANALYSIS = OUT_DIR / "analysis_result.json"    # 统计结果（相关系数等）
