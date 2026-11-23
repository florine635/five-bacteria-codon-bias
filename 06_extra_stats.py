# -*- coding: utf-8 -*-
"""
06_extra_stats.py —— 报告讨论部分用到的补充统计
==============================================

为什么单独拆一个脚本：
    报告 3.2 / 3.3 / 3.4 / 3.6 和讨论里引用了若干"派生指标"
    （相对期望曲线的残差、高表达基因与全基因池的对比、物种×聚类交叉表、
    GC3s 的组内离散、GC3s 与整体 GC 的差值……）。
    这些数不是 02/03 的默认输出，如果只存在于作者的临时笔记本里，
    读者就无法复现。所以统一放在这里，输出到 output/extra_analysis.json。

输入：output/ 下 02、03 的结果 + data/ 下的 CDS 缓存（只是为了算整体 GC）
输出：output/extra_analysis.json（同时打印到屏幕）

运行：
    python 06_extra_stats.py
"""

import json
import sys

import numpy as np
import pandas as pd

import config
import utils

# 中文 Windows 控制台是 GBK 编码，遇到无法编码的字符（如 Unicode 减号）、
# emoji 会直接抛 UnicodeEncodeError 让脚本崩掉。这里把编码错误降级为替换字符，
# 保证"打印不出来"不会变成"程序挂掉"。
try:
    sys.stdout.reconfigure(errors="replace")
except Exception:                                     # noqa: BLE001
    pass

# ---------------------------------------------------------------------------
# 高表达基因清单（HEG, highly expressed genes）
#
# 选法：按功能挑出"细菌里公认高表达"的几类——
#   核糖体蛋白（rps*/rpl*）、翻译因子（tuf/fusA/infB）、
#   糖酵解酶（gapA/pgk/eno）、分子伴侣（groL）、ATP 合酶（atpD）。
# 匹配方式是拿 NCBI [gene=...] 标签做精确匹配，命中数会因物种注释而不同，
# 这一点在报告里要写清楚（不是全覆盖，只是个方便对照的小样本）。
# ---------------------------------------------------------------------------
HEG_GENES = [
    "rpsA", "rpsB", "rpsC", "rpsI", "rplA",
    "tufA", "tufB", "fusA", "infB",
    "gapA", "gapB", "pgk", "eno", "atpD", "groL", "rpoB",
]

ORG = config.SPECIES_ORDER


def residual_analysis(df):
    """基因相对 ENC 期望曲线的残差：越负 = 选择越强。"""
    df = df.copy()
    df["enc_expected"] = utils.enc_expected(df["gc3s"].to_numpy(float))
    df["residual"] = df["enc"] - df["enc_expected"]

    out = {
        "n_genes": int(len(df)),
        "below_curve_pct": round(float((df.residual < 0).mean() * 100), 2),
        "residual_median": round(float(df.residual.median()), 2),
        "residual_min": round(float(df.residual.min()), 2),
        "residual_max": round(float(df.residual.max()), 2),
        "per_species": {},
        "strongest_selection": [],
        "closest_to_curve": [],
    }
    for org in ORG:
        g = df[df.organism == org]
        out["per_species"][org] = {
            "residual_median": round(float(g.residual.median()), 2),
            "below_curve_pct": round(float((g.residual < 0).mean() * 100), 1),
            "n": int(len(g)),
        }
    for _, r in df.nsmallest(5, "residual").iterrows():
        out["strongest_selection"].append({
            "organism": r.organism, "gene": r.gene, "n_codons": int(r.n_codons),
            "gc3s": round(float(r.gc3s), 4), "enc": round(float(r.enc), 2),
            "residual": round(float(r.residual), 2),
        })
    for _, r in df.nlargest(5, "residual").iterrows():
        out["closest_to_curve"].append({
            "organism": r.organism, "gene": r.gene, "n_codons": int(r.n_codons),
            "gc3s": round(float(r.gc3s), 4), "enc": round(float(r.enc), 2),
            "residual": round(float(r.residual), 2),
        })
    return df, out


def heg_analysis(df):
    """高表达基因 vs 全基因池：ENC 差、密码子偏好偏移。"""
    out = {"heg_list": HEG_GENES, "per_species": {}, "ecoli_codon_shift": {}}
    for org in ORG:
        g = df[df.organism == org]
        sub = g[g.gene.isin(HEG_GENES)]
        if len(sub) < 3:
            out["per_species"][org] = {"n_heg": int(len(sub)), "note": "命中太少，不比较"}
            continue
        out["per_species"][org] = {
            "n_heg": int(len(sub)),
            "enc_all": round(float(g.enc.mean()), 2),
            "enc_heg": round(float(sub.enc.mean()), 2),
            "enc_gap": round(float(g.enc.mean() - sub.enc.mean()), 2),
            "gc3s_all": round(float(g.gc3s.mean()), 4),
            "gc3s_heg": round(float(sub.gc3s.mean()), 4),
        }

    # 大肠杆菌：HEG 相对全池的 RSCU 偏移
    rcols = [c for c in df.columns if c.startswith("RSCU_")]
    ec = df[df.organism == "E. coli"]
    heg = ec[ec.gene.isin(HEG_GENES)]
    if len(heg) >= 3:
        delta = (heg[rcols].mean() - ec[rcols].mean()).sort_values(ascending=False)
        out["ecoli_codon_shift"] = {
            "n_heg": int(len(heg)),
            "up": [[c[5:], round(float(heg[c].mean()), 2),
                    round(float(ec[c].mean()), 2), round(float(v), 2)]
                   for c, v in delta.head(5).items()],
            "down": [[c[5:], round(float(heg[c].mean()), 2),
                      round(float(ec[c].mean()), 2), round(float(v), 2)]
                     for c, v in delta.tail(5).items()],
        }
    return out


def arginine_table(df):
    """精氨酸 6 个密码子的种间 RSCU（说明"同一密码子在不同物种命运相反"）。"""
    cols = ["RSCU_CGT", "RSCU_CGC", "RSCU_CGA", "RSCU_CGG", "RSCU_AGA", "RSCU_AGG"]
    m = df.groupby("organism")[cols].mean().reindex(ORG).round(2)
    return {org: {c[5:]: float(m.loc[org, c]) for c in cols} for org in ORG}


def cluster_crosstab(pca):
    """物种 × 聚类 交叉表：看错分的基因跑到了哪个簇。"""
    ct = pd.crosstab(pca.organism, pca.cluster).reindex(ORG)
    pct = (ct.div(ct.sum(axis=1), axis=0) * 100).round(1)
    dominant = ct.idxmax(axis=1)
    out = {"engine_note": "簇编号依赖 K-means 实现，以本次输出为准",
           "table_pct": {org: {f"C{c}": float(pct.loc[org, c]) for c in pct.columns}
                         for org in ORG},
           "dominant": {org: {"cluster": f"C{dominant[org]}",
                              "pct": float(pct.loc[org, dominant[org]])} for org in ORG}}
    # 每个物种最大的"误分去向"
    for org in ORG:
        others = pct.loc[org].drop(dominant[org]).sort_values(ascending=False)
        out["dominant"][org]["top_other"] = (
            {"cluster": f"C{others.index[0]}", "pct": float(others.iloc[0])}
            if len(others) else None)
    return out


def gc_dispersion(df):
    """GC3s 的组内离散，以及 GC3s 与整体 GC 的差值（偏斜）。"""
    rows = {}
    for org in ORG:
        g = df[df.organism == org].gc3s
        q1, q3 = np.percentile(g, [25, 75])
        rows[org] = {
            "mean": round(float(g.mean()), 4),
            "std": round(float(g.std()), 4),
            "cv": round(float(g.std() / g.mean()), 4),
            "iqr": round(float(q3 - q1), 4),
        }

    # 整体 CDS GC（实测，不引用外部数值）
    gc_all = {}
    for org in ORG:
        acc = config.SPECIES[org]["acc"]
        path = config.DATA_DIR / f"{acc}_cds.fna"
        if not path.exists():
            gc_all[org] = None
            continue
        gc = tot = 0
        for _, seq in utils.read_fasta(path):
            s = seq.upper()
            gc += s.count("G") + s.count("C")
            tot += len(s)
        gc_all[org] = round(100 * gc / tot, 2)

    skew = {}
    for org in ORG:
        if gc_all[org] is None:
            continue
        skew[org] = round(100 * rows[org]["mean"] - gc_all[org], 2)

    means = np.array([rows[o]["mean"] * 100 for o in ORG])
    stds = np.array([rows[o]["std"] * 100 for o in ORG])
    cvs = np.array([rows[o]["cv"] * 100 for o in ORG])
    gcv = np.array([gc_all[o] for o in ORG if gc_all[o] is not None])
    skv = np.array([skew[o] for o in ORG if o in skew])

    return {
        "per_species": rows,
        "cds_gc_pct": gc_all,
        "gc3s_minus_gc": skew,
        "corr_mean_std": round(float(np.corrcoef(means, stds)[0, 1]), 4),
        "corr_mean_cv": round(float(np.corrcoef(means, cvs)[0, 1]), 4),
        "corr_gc_skew": round(float(np.corrcoef(gcv, skv)[0, 1]), 4) if len(gcv) == len(skv) else None,
    }


def corr_robustness(df):
    """物种内 GC3s–ENC 相关，以及枯草芽孢杆菌正相关的稳健性检验。"""
    out = {"per_species": {}}
    for org in ORG:
        g = df[df.organism == org]
        out["per_species"][org] = {"r": round(float(np.corrcoef(g.gc3s, g.enc)[0, 1]), 4),
                                   "n": int(len(g))}
    b = df[df.organism == "B. subtilis"]
    out["bsubis"] = {
        "r_all": round(float(np.corrcoef(b.gc3s, b.enc)[0, 1]), 4),
        "r_gc3s_len": round(float(np.corrcoef(b.gc3s, b.n_codons)[0, 1]), 4),
        "after_min_len": {},
    }
    for thr in (150, 200, 250, 300):
        s = b[b.n_codons >= thr]
        out["bsubis"]["after_min_len"][f">={thr}"] = {
            "n": int(len(s)), "r": round(float(np.corrcoef(s.gc3s, s.enc)[0, 1]), 4)}
    return out


def misassignment_profile(df, pca):
    """
    被分到"别的物种簇"的基因是什么样？

    检验一个可证伪的推断：如果错分是因为这些基因的密码子谱介于两个物种之间，
    那么它们被分去的那一侧，应该正是 GC3s 更接近的物种。
    pca_coordinates.csv 与 codon_metrics_per_CDS.csv 行序一致（均由 03 输出），
    所以这里按位置对齐，并对长度做校验。
    """
    if len(df) != len(pca):
        return {"error": f"行数不一致 df={len(df)} pca={len(pca)}，无法对齐"}

    x = pd.DataFrame({
        "organism": df.organism.to_numpy(),
        "gene": df.gene.fillna("").to_numpy(),
        "gc3s": df.gc3s.to_numpy(float),
        "enc": df.enc.to_numpy(float),
        "cluster": pca.cluster.to_numpy(),
    })
    dominant = x.groupby("organism").cluster.agg(lambda s: s.value_counts().idxmax())
    # 反查：簇编号 -> 主导该簇的物种
    cluster_owner = {int(c): o for o, c in dominant.items()}
    species_gc3 = x.groupby("organism").gc3s.mean()

    out = {"dominant_cluster": {o: int(c) for o, c in dominant.items()}, "per_species": {}}
    for org in ORG:
        g = x[x.organism == org]
        own = g[g.cluster == dominant[org]]
        wrong = g[g.cluster != dominant[org]]
        item = {
            "n": int(len(g)), "n_misassigned": int(len(wrong)),
            "gc3s_correct": round(float(own.gc3s.mean()), 4) if len(own) else None,
            "gc3s_misassigned": round(float(wrong.gc3s.mean()), 4) if len(wrong) else None,
            "enc_correct": round(float(own.enc.mean()), 2) if len(own) else None,
            "enc_misassigned": round(float(wrong.enc.mean()), 2) if len(wrong) else None,
            "assigned_to": {},
        }
        if len(wrong):
            vc = wrong.cluster.value_counts()
            for c, n in vc.items():
                owner = cluster_owner.get(c)
                item["assigned_to"][f"C{c}"] = {
                    "n": int(n),
                    "owner_species": str(owner),
                    "owner_gc3s": round(float(species_gc3[owner]), 4) if owner in species_gc3 else None,
                    # 与接收方物种平均 GC3s 的距离，对比与自己物种的距离
                    "dist_to_owner": (round(abs(float(wrong[wrong.cluster == c].gc3s.mean())
                                               - float(species_gc3[owner])), 4)
                                      if owner in species_gc3 else None),
                    "dist_to_self": round(abs(float(wrong[wrong.cluster == c].gc3s.mean())
                                              - float(species_gc3[org])), 4),
                }
        out["per_species"][org] = item
    return out


def main():
    df = pd.read_csv(config.CSV_CDS)
    pca = pd.read_csv(config.CSV_PCA)
    print(f"读入 {len(df)} 条 CDS\n")

    df2, res = residual_analysis(df)
    print("[1] 相对 ENC 期望曲线")
    print(f"    曲线下方占比 {res['below_curve_pct']}%  "
          f"残差中位数 {res['residual_median']}  "
          f"区间 [{res['residual_min']}, {res['residual_max']}]")
    for org, v in res["per_species"].items():
        print(f"    {org:16s} 中位残差={v['residual_median']:7.2f}  下方={v['below_curve_pct']:5.1f}%")
    print("    选择最强:", [(d["gene"], d["residual"]) for d in res["strongest_selection"][:5]])

    heg = heg_analysis(df)
    print("\n[2] 高表达基因 vs 全基因池")
    for org, v in heg["per_species"].items():
        if "enc_gap" in v:
            print(f"    {org:16s} n={v['n_heg']:2d}  全池ENC={v['enc_all']:5.2f}  "
                  f"HEG={v['enc_heg']:5.2f}  差={v['enc_gap']:5.2f}")
    print("    E. coli HEG 中升高最多:", heg["ecoli_codon_shift"]["up"])

    arg = arginine_table(df)
    print("\n[3] 精氨酸密码子 RSCU")
    for org in ORG:
        print(f"    {org:16s}", arg[org])

    ct = cluster_crosstab(pca)
    print("\n[4] 物种 × 聚类")
    for org in ORG:
        print(f"    {org:16s} 主导 {ct['dominant'][org]['cluster']} "
              f"({ct['dominant'][org]['pct']}%)   "
              f"最大误分→ {ct['dominant'][org]['top_other']}")

    disp = gc_dispersion(df2)
    print("\n[5] GC3s 离散与偏斜")
    for org in ORG:
        v = disp["per_species"][org]
        print(f"    {org:16s} std={v['std']:.4f} IQR={v['iqr']:.4f} CV={v['cv']:.4f}  "
              f"CDS_GC={disp['cds_gc_pct'][org]}%  GC3s-GC={disp['gc3s_minus_gc'][org]:+.2f}")
    print(f"    r(平均GC3s, std)={disp['corr_mean_std']}  r(平均GC3s, CV)={disp['corr_mean_cv']}  "
          f"r(GC, GC3s−GC)={disp['corr_gc_skew']}")

    rob = corr_robustness(df2)
    print("\n[6] 物种内相关与稳健性")
    for org, v in rob["per_species"].items():
        print(f"    {org:16s} r={v['r']:+.4f}  n={v['n']}")
    print(f"    B. subtilis 全量 r={rob['bsubis']['r_all']:+.4f}; "
          f"r(GC3s,长度)={rob['bsubis']['r_gc3s_len']:+.4f}")
    for k, v in rob["bsubis"]["after_min_len"].items():
        print(f"      长度 {k}: n={v['n']:5d}  r={v['r']:+.4f}")

    mis = misassignment_profile(df2, pd.read_csv(config.CSV_PCA))
    print("\n[7] 错分基因的中间型检验")
    for org, v in mis.get("per_species", {}).items():
        print(f"    {org:16s} 错分 {v['n_misassigned']:4d}/{v['n']:4d}  "
              f"GC3s 正确={v['gc3s_correct']} 错分={v['gc3s_misassigned']}")
        for c, d in v["assigned_to"].items():
            if d["n"] >= 20:
                print(f"        -> {c}({d['owner_species']}, GC3s={d['owner_gc3s']}): "
                      f"{d['n']} 条  距接收方={d['dist_to_owner']}  距本物种={d['dist_to_self']}")

    result = {"residual": res, "heg": heg, "arginine_rscu": arg,
              "cluster_crosstab": ct, "gc_dispersion": disp,
              "corr_robustness": rob, "misassignment": mis}
    config.OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = config.OUT_DIR / "extra_analysis.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n已写入 {path}")


if __name__ == "__main__":
    main()
