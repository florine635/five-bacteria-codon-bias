# -*- coding: utf-8 -*-
"""
04_plot.py —— 出 4 张图
=======================

    fig1_ENC_plot.png     GC3s vs ENC 散点 + 无选择期望曲线   ← 本项目的核心图
    fig2_RSCU_heatmap.png 物种 × 61 密码子的 RSCU 热图
    fig3_PCA.png          59 维 RSCU 谱的 PCA 散点（按物种着色）
    fig4_GC3s_box.png     GC3s 的箱形图（对应课件"利用箱形图进行数据分析"）

每个 plot_* 函数都接受外部传入的 ax，所以 05_gui.py 可以直接复用来画图，
不用把代码抄两遍。
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import config
import utils

utils.setup_chinese_font()

ORG = config.SPECIES_ORDER
COL = config.COLORS


# ---------------------------------------------------------------- 图1
def plot_enc(ax, df, orgs=None, n_per_org=400, seed=42):
    """ENC-plot：横轴 GC3s，纵轴 ENC，虚线是"仅突变压力"的期望曲线。"""
    orgs = orgs or ORG
    xs = np.linspace(0.005, 0.995, 400)
    ax.plot(xs, utils.enc_expected(xs), "k--", lw=1.3,
            label="期望曲线（无选择）", zorder=1)

    rng = np.random.default_rng(seed)
    for org in orgs:
        g = df[df.organism == org]
        if g.empty:
            continue
        idx = rng.choice(len(g), size=min(n_per_org, len(g)), replace=False)
        ax.scatter(g.gc3s.iloc[idx], g.enc.iloc[idx], s=11, alpha=.55,
                   c=COL.get(org, "gray"), edgecolors="none", zorder=2,
                   label=f"{org} (n={len(g)})")

    ax.set_xlabel("GC3s　第三位密码子的 GC 含量")
    ax.set_ylabel("ENC　有效密码子数")
    ax.set_title("ENC-plot：5 个细菌基因组的密码子使用偏好", fontsize=11)
    ax.set_xlim(0, 1)
    ax.set_ylim(20, 62)
    ax.grid(alpha=.25, ls=":")
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")


# ---------------------------------------------------------------- 图2
def plot_rscu_heatmap(fig, ax, rscu_matrix, table=None):
    """
    RSCU 热图：行=物种，列=密码子（按氨基酸分组排列，组间画分隔线）。

    rscu_matrix: index 是物种，columns 是 RSCU_XXX
    """
    if table is None:
        table = utils.get_table(11)

    codons = [c[5:] for c in rscu_matrix.columns]
    order = sorted(codons, key=lambda c: (table.forward_table.get(c, "Z"), c))
    mat = rscu_matrix[["RSCU_" + c for c in order]].to_numpy(float)

    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0, vmax=2)
    ax.set_yticks(range(len(rscu_matrix)))
    ax.set_yticklabels(list(rscu_matrix.index), fontsize=8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=5.5)

    # 氨基酸分组分隔线
    aas = [table.forward_table.get(c, "Z") for c in order]
    for i in range(1, len(aas)):
        if aas[i] != aas[i - 1]:
            ax.axvline(i - 0.5, color="white", lw=.7, alpha=.8)

    ax.set_title("相对同义密码子使用度 RSCU（蓝 < 1 回避，红 > 1 偏好）",
                 fontsize=9.5)
    cb = fig.colorbar(im, ax=ax, pad=.008, fraction=.015)
    cb.ax.tick_params(labelsize=6)
    return im


# ---------------------------------------------------------------- 图3
def plot_pca(ax, pca_df, orgs=None, n_max=2500, seed=7, show_cluster=False,
             explained=None):
    """
    PCA 散点：按物种着色。

    explained: 可选，前两个主成分的解释方差百分比（如 [39.89, 8.41]）。
               有就写进轴标签 —— 报告图里这是必要信息，缺了审稿会问。
    """
    orgs = orgs or ORG
    rng = np.random.default_rng(seed)
    if len(pca_df) > n_max:
        idx = rng.choice(len(pca_df), size=n_max, replace=False)
        sub = pca_df.iloc[idx]
    else:
        sub = pca_df

    for org in orgs:
        m = sub.organism == org
        if not m.any():
            continue
        ax.scatter(sub.PC1[m], sub.PC2[m], s=9, alpha=.6, c=COL.get(org, "gray"),
                   edgecolors="none", label=org)
    if show_cluster and "cluster" in sub.columns:
        for cid in sorted(sub.cluster.unique()):
            m = sub.cluster == cid
            ax.annotate(f"C{cid}", (sub.PC1[m].mean(), sub.PC2[m].mean()),
                        fontsize=8, ha="center", va="center", color="black",
                        bbox=dict(boxstyle="round,pad=.25", fc="white",
                                  ec="gray", alpha=.85))
    ax.set_xlabel(f"PC1（{explained[0]:.1f}%）" if explained else "PC1")
    ax.set_ylabel(f"PC2（{explained[1]:.1f}%）" if explained else "PC2")
    ax.set_title("RSCU 谱的主成分分析", fontsize=10.5)
    ax.legend(fontsize=7.5, frameon=False, markerscale=2)
    ax.grid(alpha=.2, ls=":")


# ---------------------------------------------------------------- 图4
def plot_box(ax, df, orgs=None):
    """GC3s 箱形图：一眼看出组间差异与组内离散。"""
    orgs = orgs or ORG
    data, labels, colors = [], [], []
    for org in orgs:
        g = df[df.organism == org]
        if g.empty:
            continue
        data.append(g.gc3s.to_numpy())
        labels.append(org)
        colors.append(COL.get(org, "gray"))

    bp = ax.boxplot(data, patch_artist=True, widths=.55, showfliers=False,
                    medianprops=dict(color="black", lw=1.2))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(.75)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=15, fontsize=8)
    ax.set_ylabel("GC3s")
    ax.set_title("各物种 GC3s 分布（箱线图）", fontsize=10.5)
    ax.grid(axis="y", alpha=.25, ls=":")


def load_explained():
    """从 03 的输出里读 PCA 解释方差百分比；读不到就返回 None（轴标签自动降级）。"""
    import json
    try:
        j = json.loads(config.JSON_ANALYSIS.read_text(encoding="utf-8"))
        return j["pca"]["explained_pct"][:2]
    except Exception:                                  # noqa: BLE001
        return None


# ---------------------------------------------------------------- 主流程
def main():
    df = pd.read_csv(config.CSV_CDS)
    rscu_matrix = pd.read_csv(config.CSV_RSCU, index_col=0)
    pca_df = pd.read_csv(config.CSV_PCA)
    explained = load_explained()
    print(f"读入 {len(df)} 条 CDS")

    # 图1
    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=220)
    plot_enc(ax, df)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig1_ENC_plot.png")
    plt.close(fig)

    # 图2
    fig, ax = plt.subplots(figsize=(12.5, 2.6), dpi=220)
    plot_rscu_heatmap(fig, ax, rscu_matrix)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig2_RSCU_heatmap.png")
    plt.close(fig)

    # 图3
    fig, ax = plt.subplots(figsize=(6.4, 5.2), dpi=220)
    plot_pca(ax, pca_df, show_cluster=True, explained=explained)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig3_PCA.png")
    plt.close(fig)

    # 图4
    fig, ax = plt.subplots(figsize=(6.6, 4.8), dpi=220)
    plot_box(ax, df)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig4_GC3s_box.png")
    plt.close(fig)

    print(f"4 张图已写入 {config.FIG_DIR}")
    print("下一步 → python 05_gui.py（图形界面）")


if __name__ == "__main__":
    main()
