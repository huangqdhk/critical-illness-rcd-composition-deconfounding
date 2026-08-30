# -*- coding: utf-8 -*-
"""
run_demo_all.py — 一键运行两个采纳演示（发布版入口）
========================================================================
用法：python run_demo_all.py [--outdir 目录]
两个演示数据集的输入矩阵由用户提供（GEO 官网下载，几 MB 级 processed 文件），
路径通过 --gse157103 / --gse66099 传入；缺省时使用项目本地路径。
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "demos"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import demo_GSE157103 as d1  # noqa: E402
import demo_GSE66099 as d2  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse157103", default=None, help="GSE157103_genes.tpm.tsv.gz 路径")
    ap.add_argument("--gse66099", default=None, help="GSE66099_series_matrix.txt.gz 路径")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    print("== 采纳演示 1：GSE157103", flush=True)
    d1.main(a.gse157103, a.outdir)
    print("\n== 采纳演示 2：GSE66099", flush=True)
    d2.main(a.gse66099, a.outdir)
    print("\n== 全部完成", flush=True)


if __name__ == "__main__":
    main()
