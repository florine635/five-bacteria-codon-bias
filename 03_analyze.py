# -*- coding: utf-8 -*-
"""
03_analyze.py —— 统计与降维聚类
==============================

做了什么：
    1. GC3s 与 ENC 的相关性（全基因池、以及每个物种内部各算一次）
    2. 对 59 维 RSCU 谱做主成分分析（PCA，用 SVD 手写，不依赖 sklearn）
    3. K-means 聚类，并算"聚类纯度"——看聚出来的类是否等于物种
    4. 结果写 output/analysis_result.json 与 output/pca_coordinates.csv

注意：PCA 是手写 SVD，sklearn 只用于 K-means；万一没装 sklearn，
会自动退回到自带的 numpy 版 K-means，脚本不会中断。
"""

import json

import numpy as np
import pandas as pd

import config


# ---------------------------------------------------------------- 工具
def pearson(x, y):
    """皮尔逊相关系数（空数组返回 nan）。"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan")
    x, y = x - x.mean(), y - y.mean()
    denom = np.sqrt((x ** 2).sum() * (y ** 2).sum())
    return float((x * y).sum() / denom) if denom else float("nan")


def pca_svd(X, k=2):
    """
    PCA（手写版）：标准化 → SVD → 取前 k 个主成分。

    为什么不用 sklearn.decomposition.PCA？
    自己写一遍才说得清"主成分就是协方差矩阵的特征向量"，
    而且这里只要 2 维，SVD 足够。
    """
    Xc = X - X.mean(axis=0)
    # 标准差为 0 的列（该密码子在所有基因里都一样）避免除零
    sd = Xc.std(axis=0)
    sd[sd == 0] = 1.0
    Xz = Xc / sd

    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    scores = U[:, :k] * S[:k]
    explained = (S ** 2) / (S ** 2).sum()          # 每个主成分的解释方差比例
    # 前两个主成分在各密码子上的载荷，用于解释"PC1 代表什么"
    loadings = Vt[:k]
    return scores, explained, loadings


def kmeans_fallback(X, k=5, n_iter=300, n_init=10, seed=42):
    """
    numpy 版 K-means（sklearn 不可用时的后备）。返回标签数组。

    为什么要写这么长：朴素随机初始化 + 单次运行会让结果**依赖有没有装 sklearn**，
    同一份数据可能得到 0.9143 或 0.9145 两种纯度，报告数字就不稳了。
    所以这里照 sklearn 的做法实现了 k-means++ 初始化 + 多次重启取最优，
    使两种引擎收敛到同一个解。
    """
    rng = np.random.default_rng(seed)
    best_labels, best_inertia = None, np.inf

    for _ in range(n_init):
        # ---- k-means++ 初始化：让初始中心彼此尽量远 ----
        centers = [X[rng.integers(len(X))]]
        d2 = ((X - centers[0]) ** 2).sum(axis=1)
        for _k in range(1, k):
            total = d2.sum()
            probs = np.full(len(X), 1.0 / len(X)) if total <= 0 else d2 / total
            idx = int(rng.choice(len(X), p=probs))
            centers.append(X[idx])
            d2 = np.minimum(d2, ((X - X[idx]) ** 2).sum(axis=1))
        centers = np.array(centers)

        # ---- Lloyd 迭代 ----
        labels = np.full(len(X), -1, dtype=int)
        for _it in range(n_iter):
            d = ((X[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
            new = d.argmin(axis=1)
            if np.array_equal(new, labels):
                break
            labels = new
            for j in range(k):
                if (labels == j).any():
                    centers[j] = X[labels == j].mean(axis=0)

        inertia = float(((X - centers[labels]) ** 2).sum())
        if inertia < best_inertia:
            best_inertia, best_labels = inertia, labels

    return best_labels


def purity_score(true_labels, cluster_labels):
    """
    聚类纯度：每个簇取"占多数的真实类别"，正确数 / 总数。

    纯度 = 1 表示聚类结果与物种完全一致；但纯度会随簇数增加而虚高，
    所以报告里要同时给出簇数（这里固定 5，正好等于物种数）。
    """
    total, correct = 0, 0
    detail = []
    for j in np.unique(cluster_labels):
        mask = cluster_labels == j
        vals, counts = np.unique(true_labels[mask], return_counts=True)
        top = counts.argmax()
        correct += counts[top]
        total += mask.sum()
        detail.append({
            "cluster": int(j),
            "size": int(mask.sum()),
            "majority_species": str(vals[top]),
            "majority_pct": round(100 * counts[top] / mask.sum(), 2),
            "composition": {str(v): int(c) for v, c in zip(vals, counts)},
        })
    return correct / total if total else float("nan"), detail


# ---------------------------------------------------------------- 主流程
def main():
    df = pd.read_csv(config.CSV_CDS)
    rscu_cols = [c for c in df.columns if c.startswith("RSCU_")]
    print(f"读入 {len(df)} 条 CDS，{len(rscu_cols)} 个 RSCU 维度")

    result = {}

    # ---- 1. 相关性 ----
    r_all = pearson(df["gc3s"], df["enc"])
    result["corr_gc3s_enc_all"] = {"r": round(r_all, 4), "n": int(len(df))}
    print(f"\n[1] GC3s 与 ENC 相关（全部 {len(df)} 条）: r = {r_all:.4f}")

    per_org = {}
    for name in config.SPECIES_ORDER:
        g = df[df.organism == name]
        per_org[name] = {"r": round(pearson(g["gc3s"], g["enc"]), 4),
                         "n": int(len(g)),
                         "gc3s_mean": round(float(g.gc3s.mean()), 4),
                         "enc_mean": round(float(g.enc.mean()), 2)}
        print(f"    {name:<17} r = {per_org[name]['r']:>7.4f}  n = {per_org[name]['n']}")
    result["corr_per_species"] = per_org

    # 物种层面（每物种一个点，n=5）——样本量小，只作描述
    sm = df.groupby("organism").agg(gc3s=("gc3s", "mean"), enc=("enc", "mean"))
    sm = sm.reindex(config.SPECIES_ORDER)
    r_species = pearson(sm["gc3s"], sm["enc"])
    result["corr_species_level"] = {"r": round(r_species, 4), "n": int(len(sm))}
    print(f"[1b] 物种层面（n={len(sm)}）r = {r_species:.4f}  ← 样本太少，仅描述")

    # ---- 2. PCA ----
    X = df[rscu_cols].to_numpy(float)
    scores, explained, loadings = pca_svd(X, k=2)
    result["pca"] = {
        "explained_pct": [round(float(e) * 100, 2) for e in explained[:5]],
        "n_features": len(rscu_cols),
    }
    print(f"\n[2] PCA 前两个主成分解释方差: "
          f"PC1 = {explained[0] * 100:.2f}%, PC2 = {explained[1] * 100:.2f}%")

    # PC1 载荷最大的密码子（帮助解释 PC1 代表什么）
    pc1 = pd.Series(loadings[0], index=[c[5:] for c in rscu_cols])
    result["pca"]["pc1_top_positive"] = [[c, round(float(v), 3)]
                                         for c, v in pc1.nlargest(5).items()]
    result["pca"]["pc1_top_negative"] = [[c, round(float(v), 3)]
                                         for c, v in pc1.nsmallest(5).items()]
    print("    PC1 正端载荷最高密码子:", result["pca"]["pc1_top_positive"])
    print("    PC1 负端载荷最高密码子:", result["pca"]["pc1_top_negative"])

    pca_df = pd.DataFrame({
        "organism": df["organism"],
        "gene": df["gene"],
        "PC1": scores[:, 0],
        "PC2": scores[:, 1],
    })
    pca_df.to_csv(config.CSV_PCA, index=False, encoding="utf-8-sig")

    # ---- 3. K-means ----
    Xz = (X - X.mean(0)) / np.where(X.std(0) == 0, 1.0, X.std(0))
    used = "sklearn"
    try:
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=len(config.SPECIES_ORDER), n_init=10, random_state=42)
        labels = km.fit_predict(Xz)
    except Exception as exc:                          # noqa: BLE001
        used = f"numpy-fallback ({type(exc).__name__})"
        labels = kmeans_fallback(Xz, k=len(config.SPECIES_ORDER))
    pur, detail = purity_score(df["organism"].to_numpy(), labels)
    result["kmeans"] = {
        "engine": used,
        "k": len(config.SPECIES_ORDER),
        "purity": round(float(pur), 4),
        "clusters": detail,
    }
    print(f"\n[3] K-means（{used}，k={len(config.SPECIES_ORDER)}）聚类纯度 = {pur:.4f}")
    for d in detail:
        print(f"    簇 {d['cluster']}: {d['size']:>5} 条，"
              f"主导 {d['majority_species']} ({d['majority_pct']}%)")

    pca_df["cluster"] = labels
    pca_df.to_csv(config.CSV_PCA, index=False, encoding="utf-8-sig")

    config.JSON_ANALYSIS.write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n结果已写入 {config.JSON_ANALYSIS}")
    print(f"[注意] 报告里请用**本次**输出的簇编号与条数（引擎：{used}），"
          f"不同 K-means 实现会让个别边界基因换簇。")
    print("下一步 → python 04_plot.py")


if __name__ == "__main__":
    main()
