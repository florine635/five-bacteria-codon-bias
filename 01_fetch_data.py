# -*- coding: utf-8 -*-
"""
01_fetch_data.py —— 从 NCBI 取数并缓存
=====================================

用法：
    python 01_fetch_data.py --verify     # 只核实登录号对不对（强烈建议先跑一次）
    python 01_fetch_data.py              # 下载全部 CDS 到 data/

做了什么：
    1. 对 config.SPECIES 里每个染色体登录号，向 NCBI 查摘要，打印物种名和长度，
       让你**用眼睛确认**没有填错登录号（这一步是防翻车的关键）。
    2. 用 rettype="fasta_cds_na" 一次性取回该基因组注释过的全部 CDS
       （比逐条取快得多，而且拿到的是"全部 CDS"，不是随手挑的样本）。
    3. 结果缓存到 data/<登录号>_cds.fna。下次运行直接读缓存，
       所以**演示时断网也能跑**。

为什么用染色体登录号 NC_ 而不是 GCF_ 组装号：
    NCBI 的 nuccore 库不认 GCF_（会报 "Failed to understand id"），
    必须用染色体自身的 NC_ 登录号。这是实测踩过的坑。
"""

import argparse
import sys
import time

try:
    from Bio import Entrez
except ImportError:
    print("这个脚本需要 Biopython —— 它用 Bio.Entrez 去 NCBI 取数。\n"
          "\n安装（清华镜像，很快）：\n"
          "    pip install -i https://pypi.tuna.tsinghua.edu.cn/simple biopython\n"
          "\n提示：如果你只是想在**已缓存的数据**上算指标、画图或打开界面，\n"
          "      02 / 03 / 04 / 05 都不需要 Biopython，可以直接运行。\n")
    sys.exit(1)

import config
import utils


def _retry(fn, tries=None, label=""):
    """带重试地调用一个函数；NCBI 偶发连接重置，重试基本都能过。"""
    tries = tries or config.NCBI_RETRY
    last = None
    for i in range(tries):
        try:
            return fn()
        except Exception as exc:                      # noqa: BLE001
            last = exc
            print(f"    [重试 {i + 1}/{tries}] {label} {type(exc).__name__}")
            time.sleep(2 + i)
    raise RuntimeError(f"{label} 连续失败 {tries} 次：{last}")


def verify_accessions():
    """核实登录号：打印 NCBI 记录的物种名与长度。"""
    Entrez.email = config.NCBI_EMAIL
    accs = [v["acc"] for v in config.SPECIES.values()]

    def _fetch():
        handle = Entrez.esummary(db="nuccore", id=",".join(accs))
        docs = Entrez.read(handle)
        handle.close()
        return docs

    docs = _retry(_fetch, label="esummary")
    print(f"{'登录号':<16}{'长度(bp)':>12}  标题")
    print("-" * 90)
    for d in docs:
        acc = str(d.get("AccessionVersion", "?"))
        print(f"{acc:<16}{int(d.get('Length', 0)):>12,}  {str(d.get('Title', ''))[:60]}")
    print("-" * 90)
    print("请逐行核对：物种名、长度是否与你预期一致。不一致就改 config.SPECIES。")


def fetch_one(name, acc):
    """取一个基因组的全部 CDS，优先用缓存。返回 (记录数, 来源)。"""
    cache = config.DATA_DIR / f"{acc}_cds.fna"
    if cache.exists() and cache.stat().st_size > 5000:
        n = sum(1 for line in open(cache, encoding="utf-8", errors="replace")
                if line.startswith(">"))
        return n, "缓存"

    def _fetch():
        handle = Entrez.efetch(db="nuccore", id=acc,
                               rettype="fasta_cds_na", retmode="text")
        text = handle.read()
        handle.close()
        return text

    text = _retry(_fetch, label=f"efetch {acc}")
    time.sleep(config.NCBI_SLEEP)

    if not text.startswith(">"):
        raise RuntimeError(f"{acc} 返回内容异常：{text[:120]!r}\n"
                           f"常见原因：登录号不是 NC_ 开头的染色体号。")
    cache.write_text(text, encoding="utf-8")
    n = sum(1 for line in text.splitlines() if line.startswith(">"))
    return n, "下载"


def main():
    ap = argparse.ArgumentParser(description="下载细菌基因组 CDS")
    ap.add_argument("--verify", action="store_true", help="只核实登录号，不下载")
    args = ap.parse_args()

    if args.verify:
        verify_accessions()
        return

    Entrez.email = config.NCBI_EMAIL
    if "example.edu" in config.NCBI_EMAIL:
        print("[注意] config.NCBI_EMAIL 还是示例值，建议改成真实邮箱。\n")

    print(f"开始取数，共 {len(config.SPECIES)} 个基因组 → {config.DATA_DIR}\n")
    total = 0
    for name in config.SPECIES_ORDER:
        acc = config.SPECIES[name]["acc"]
        n, src = fetch_one(name, acc)
        total += n
        print(f"  {name:<17} {acc:<14} {n:>5} 条 CDS   [{src}]")
        time.sleep(config.NCBI_SLEEP)

    print(f"\n完成：{total} 条 CDS 已缓存。下一步 → python 02_compute_metrics.py")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:                            # noqa: BLE001
        print(f"\n[失败] {e}", file=sys.stderr)
        print("排查顺序：① 网络 ② 邮箱是否填了 ③ 登录号是否 NC_ 开头 ④ 稍后重试",
              file=sys.stderr)
        sys.exit(1)
