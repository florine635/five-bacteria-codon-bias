# -*- coding: utf-8 -*-
"""
02_compute_metrics.py —— 清洗 CDS 并计算三个指标
================================================

输入：data/<登录号>_cds.fna（01 下载的缓存）
输出：
    output/codon_metrics_per_CDS.csv   每条 CDS 一行（物种、基因、长度、GC3s、ENC、59 个 RSCU）
    output/species_summary.csv         每个物种一行（CDS 数、平均 GC3s、平均 ENC ...）
    output/rscu_matrix.csv             物种 × 密码子 的 RSCU 均值（画热图用）
    output/filter_log.csv              过滤统计：每种原因滤掉多少条（报告"方法"要用）

这一步不联网，纯计算，可以随便重跑。
"""

import json
from collections import Counter

import pandas as pd

import config
import utils


def main():
    print("开始计算指标...")
    per_cds = []
    filter_rows = []

    for name in config.SPECIES_ORDER:
        acc = config.SPECIES[name]["acc"]
        table = utils.get_table(config.SPECIES[name]["table"])
        families, _stops = utils.codon_families(table)

        path = config.DATA_DIR / f"{acc}_cds.fna"
        if not path.exists():
            raise FileNotFoundError(
                f"缺少缓存 {path}，请先运行：python 01_fetch_data.py")

        reasons = Counter()
        n_valid = 0
        n_records = 0
        for header, seq in utils.read_fasta(path):
            n_records += 1
            codons, reason = utils.clean_cds(
                seq, table,
                min_codons=config.MIN_CODONS,
                require_atg=config.REQUIRE_ATG,
                min_acgt=config.MIN_FRACTION_ACGT,
            )
            if codons is None:
                reasons[reason] += 1          # 只统计"被剔除"的原因，不统计 ok
                continue

            counts = utils.count_codons(codons)
            enc = utils.enc_value(counts, families)

            # ENC 溢出保护：超出 [20, 61] 说明该基因在某个高简并家族里密码子太少，
            # 估计不可靠 —— 剔除并记录，而不是硬截断。
            if not (utils.ENC_MIN <= enc <= utils.ENC_MAX):
                reasons["enc_out_of_range"] += 1
                continue

            n_valid += 1
            rscu = utils.rscu_vector(counts, families)
            row = {
                "organism": name,
                "genome": acc,
                "protein_id": header.split()[0] if header else "",
                "gene": utils.parse_gene(header),
                "len_bp": len(seq),
                "n_codons": len(codons),
                "gc3s": utils.gc3s(codons, table),
                "enc": enc,
            }
            for c, v in rscu.items():
                row["RSCU_" + c] = v
            per_cds.append(row)

        total = n_records
        filter_rows.append({
            "organism": name,
            "total_records": total,
            "kept": n_valid,
            "kept_pct": round(100 * n_valid / total, 2) if total else 0,
            **{f"drop_{k}": v for k, v in reasons.items()},
        })
        print(f"  {name:<17} 记录 {total:>5}  保留 {n_valid:>5}  "
              f"({100 * n_valid / total:.1f}%)" if total else f"  {name}: 无记录")

    df = pd.DataFrame(per_cds)
    if df.empty:
        raise SystemExit("没有可用的 CDS，检查 01 是否成功、过滤阈值是否过严。")

    # 补齐缺失的 RSCU 列（某物种可能没有某个密码子）
    rscu_cols = sorted([c for c in df.columns if c.startswith("RSCU_")],
                       key=lambda c: c[5:])
    df[rscu_cols] = df[rscu_cols].fillna(0.0)

    df.to_csv(config.CSV_CDS, index=False, encoding="utf-8-sig")

    # ---- 物种汇总 ----
    summary = (df.groupby("organism")
                 .agg(n_CDS=("enc", "size"),
                      n_codons_mean=("n_codons", "mean"),
                      gc3s_mean=("gc3s", "mean"),
                      gc3s_std=("gc3s", "std"),
                      enc_mean=("enc", "mean"),
                      enc_std=("enc", "std"),
                      enc_min=("enc", "min"),
                      enc_max=("enc", "max"))
                 .reindex(config.SPECIES_ORDER)
                 .round(4))
    summary.to_csv(config.CSV_SUMMARY, encoding="utf-8-sig")

    # ---- RSCU 均值矩阵（热图用）----
    rscu_matrix = df.groupby("organism")[rscu_cols].mean().reindex(config.SPECIES_ORDER)
    rscu_matrix.to_csv(config.CSV_RSCU, encoding="utf-8-sig")

    # ---- 过滤日志 ----
    # 不同物种触发的剔除原因不同，DataFrame 会有空值：填 0 并转成整数，
    # 否则 CSV 里会出现 "2.0"、空单元格这种不好引用的样子。
    flog = pd.DataFrame(filter_rows)
    for col in flog.columns:
        if col.startswith("drop_"):
            flog[col] = flog[col].fillna(0).astype(int)
    flog.to_csv(config.CSV_FILTERLOG, index=False, encoding="utf-8-sig")

    print("\n物种汇总：")
    print(summary[["n_CDS", "gc3s_mean", "enc_mean"]].to_string())
    print(f"\n输出：\n  {config.CSV_CDS}\n  {config.CSV_SUMMARY}"
          f"\n  {config.CSV_RSCU}\n  {config.CSV_FILTERLOG}")
    print("\n下一步 → python 03_analyze.py")


if __name__ == "__main__":
    main()
