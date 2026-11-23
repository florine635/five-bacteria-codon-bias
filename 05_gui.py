# -*- coding: utf-8 -*-
"""
05_gui.py —— 图形界面（tkinter，对应课件第 7 周"面向对象和图形用户界面"）
=====================================================================

运行：
    python 05_gui.py

界面能力：
    · 勾选要显示/比较的物种
    · 下拉切换 4 种图（ENC-plot / RSCU 热图 / PCA / GC3s 箱形图）
    · 「数据表」页签看每个物种的汇总指标（ttk.Treeview）
    · 可以把当前图导出成 PNG

设计说明（答辩可以讲的三点）：
    1. 界面类 CodonApp 继承 tk.Tk，把"数据"与"视图"分开：
       数据只从 output/*.csv 读一次，视图（_redraw）只负责重画。
    2. 绘图函数复用 04_plot.py，没有把代码抄一遍 —— 这就是模块化的好处。
    3. 界面**不联网、不重算**，启动即用；重活都在 02/03 里跑完了。
       这样演示时不会因为网络或计算卡住。
"""

import sys
from importlib import import_module


def _preflight():
    """
    启动前检查依赖，缺了就打印"能直接复制去执行"的修复命令，而不是甩一屏 traceback。

    注意：这里**不检查 biopython**。看图只需要密码子对照表，
    utils.py 里带了内置表，没装 Biopython 也能正常出图。
    只有 01_fetch_data.py（联网取数）才真正需要 Biopython。
    """
    checks = [
        ("tkinter", "Python 自带；Windows 版 Anaconda 一般都有，缺了就重装 Python 并勾选 tcl/tk"),
        ("pandas", "pip install pandas"),
        ("matplotlib", "pip install matplotlib"),
    ]
    missing = []
    for mod, hint in checks:
        try:
            __import__(mod)
        except ImportError:
            missing.append((mod, hint))

    if missing:
        print("缺少依赖，无法启动界面：")
        for mod, hint in missing:
            print(f"  · {mod:12s} → {hint}")
        print("\n一次性装齐（清华镜像）：")
        print("  pip install -i https://pypi.tuna.tsinghua.edu.cn/simple "
              "pandas matplotlib biopython scikit-learn")
        sys.exit(1)


_preflight()

import tkinter as tk
from tkinter import ttk, messagebox

import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import config
import utils

# 04_plot 的模块名以数字开头，不能直接 import，这里用 importlib
plot = import_module("04_plot")


class _SizeEvent:
    """给 canvas.resize() 用的最小事件对象（那份代码只读 width / height）。"""

    def __init__(self, width, height):
        self.width = width
        self.height = height


class CodonApp(tk.Tk):
    """密码子偏好分析结果浏览器。"""

    CHARTS = ["ENC-plot", "RSCU 热图", "PCA 散点", "GC3s 箱形图"]

    def __init__(self):
        super().__init__()
        self.title("细菌密码子使用偏好分析 — 结果浏览")
        self.geometry("1180x760")
        self.minsize(900, 600)

        self._load_data()
        self._build_widgets()
        self._fill_table()
        self._redraw()

    # ------------------------------------------------------------ 数据
    def _load_data(self):
        """只读 output/ 里已算好的结果，不重算、不联网。"""
        missing = [p.name for p in (config.CSV_CDS, config.CSV_RSCU,
                                    config.CSV_SUMMARY, config.CSV_PCA)
                   if not p.exists()]
        if missing:
            messagebox.showerror(
                "缺少结果文件",
                "找不到：" + "、".join(missing) +
                "\n\n请先依次运行：\n"
                "  python 01_fetch_data.py\n"
                "  python 02_compute_metrics.py\n"
                "  python 03_analyze.py\n"
                "  python 04_plot.py")
            raise SystemExit(1)

        self.df = pd.read_csv(config.CSV_CDS)
        self.rscu = pd.read_csv(config.CSV_RSCU, index_col=0)
        self.summary = pd.read_csv(config.CSV_SUMMARY, index_col=0)
        self.pca = pd.read_csv(config.CSV_PCA)
        self.orgs = [o for o in config.SPECIES_ORDER if o in set(self.df.organism)]

    # ------------------------------------------------------------ 界面
    def _build_widgets(self):
        utils.setup_chinese_font()

        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("TLabel", font=("Microsoft YaHei UI", 10))
        style.configure("TButton", font=("Microsoft YaHei UI", 10))
        style.configure("TCheckbutton", font=("Microsoft YaHei UI", 9))
        style.configure("Treeview", font=("Microsoft YaHei UI", 9), rowheight=24)
        style.configure("Treeview.Heading", font=("Microsoft YaHei UI", 9, "bold"))

        # ---- 顶部：标题 + 控制条 ----
        top = ttk.Frame(self, padding=(12, 10, 12, 6))
        top.pack(side="top", fill="x")

        ttk.Label(top, text="细菌密码子使用偏好分析",
                  font=("Microsoft YaHei UI", 15, "bold")).pack(side="left")
        ttk.Label(top, text=f"　共 {len(self.df)} 条 CDS / {len(self.orgs)} 个物种",
                  foreground="#666").pack(side="left")

        ctrl = ttk.Frame(self, padding=(12, 0, 12, 8))
        ctrl.pack(side="top", fill="x")

        ttk.Label(ctrl, text="物种：").pack(side="left")
        self.species_vars = {}
        for org in self.orgs:
            var = tk.BooleanVar(value=True)
            self.species_vars[org] = var
            ttk.Checkbutton(ctrl, text=org, variable=var,
                            command=self._redraw).pack(side="left", padx=(0, 8))

        ttk.Label(ctrl, text="　图表：").pack(side="left")
        self.chart_var = tk.StringVar(value=self.CHARTS[0])
        combo = ttk.Combobox(ctrl, textvariable=self.chart_var, values=self.CHARTS,
                             state="readonly", width=14)
        combo.pack(side="left")
        combo.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        ttk.Button(ctrl, text="全选", command=self._select_all).pack(side="left", padx=6)
        ttk.Button(ctrl, text="导出当前图", command=self._export).pack(side="left")

        # ---- 中部：页签（图表 / 数据表）----
        nb = ttk.Notebook(self)
        nb.pack(side="top", fill="both", expand=True, padx=12, pady=(0, 6))

        plot_tab = ttk.Frame(nb)
        nb.add(plot_tab, text="  图表  ")

        self.fig = Figure(figsize=(9.0, 5.6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_tab)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # 记住上一次的画布尺寸：窗口缩放时按需重画，尺寸没变就不重画
        # （否则 <Configure> 会和 _redraw 互相触发，变成死循环）
        self._last_canvas_size = (0, 0)
        self.canvas.get_tk_widget().bind("<Configure>", self._on_canvas_resize)

        table_tab = ttk.Frame(nb)
        nb.add(table_tab, text="  数据表  ")
        cols = ("n_CDS", "gc3s_mean", "enc_mean", "enc_min", "enc_max")
        self.tree = ttk.Treeview(table_tab, columns=cols, show="tree headings",
                                 height=10)
        self.tree.heading("#0", text="物种")
        self.tree.heading("n_CDS", text="CDS 数")
        self.tree.heading("gc3s_mean", text="平均 GC3s")
        self.tree.heading("enc_mean", text="平均 ENC")
        self.tree.heading("enc_min", text="ENC 最小")
        self.tree.heading("enc_max", text="ENC 最大")
        self.tree.column("#0", width=200)
        for c in cols:
            self.tree.column(c, width=110, anchor="center")
        self.tree.pack(side="left", fill="both", expand=True, padx=(0, 6), pady=8)

        note = ttk.Label(
            table_tab, justify="left", foreground="#444",
            text=("读法：\n"
                  "· GC3s 越高 → 突变压力越偏向 G/C\n"
                  "· ENC 越低 → 密码子使用越偏（20=极端偏好，61=完全随机）\n"
                  "· 两者通常负相关：GC3s 高的物种，择压也更集中\n\n"
                  "箱线图与热图在 04_plot.py 里保存了 220 dpi 版本，\n"
                  "报告里建议直接用 figures/ 下的图，清晰度更好。"))
        note.pack(side="left", fill="y", padx=8, pady=8)

        # ---- 底部状态栏 ----
        self.status = ttk.Label(self, text="就绪", anchor="w",
                                padding=(12, 4), relief="sunken")
        self.status.pack(side="bottom", fill="x")

    def _fill_table(self):
        for org in self.orgs:
            if org not in self.summary.index:
                continue
            r = self.summary.loc[org]
            self.tree.insert("", "end", text=org, values=(
                int(r["n_CDS"]),
                f"{r['gc3s_mean']:.3f}",
                f"{r['enc_mean']:.2f}",
                f"{r['enc_min']:.1f}",
                f"{r['enc_max']:.1f}",
            ))

    # ------------------------------------------------------------ 交互
    def _selected_species(self):
        sel = [o for o in self.orgs if self.species_vars[o].get()]
        return sel or self.orgs          # 全不选时退回全选，避免空白图

    def _select_all(self):
        for var in self.species_vars.values():
            var.set(True)
        self._redraw()

    def _sync_figure_size(self):
        """
        让「figure 尺寸」和「Tk 上的位图尺寸」都精确等于画布控件的尺寸。

        为什么要同时管两个（这里踩过两次坑）：
          · figure 决定渲染分辨率。若 figure 比控件小，matplotlib 只重绘
            figure 覆盖的那块像素，控件四周残留上一张图 → "多张图叠在一起"。
          · Tk 上真正显示的是 FigureCanvasTkAgg 内部那张 PhotoImage（_tkphoto）。
            它只在画布 resize 时重建；若它比 figure 小，就只有左上角一块被更新，
            右边和下边直接被裁掉 → "每张图右边少一块"。
          · 而 matplotlib 自己的 resize 处理器挂在同一个 <Configure> 事件上，
            我们一旦用 bind("<Configure>", ...) 就会把它顶掉（Tk 的 bind 默认是
            替换而非追加）。所以这里**主动调用 canvas.resize()**，
            不依赖事件绑定的先后顺序。
        """
        tk = self.canvas.get_tk_widget()
        w, h = tk.winfo_width(), tk.winfo_height()
        if w < 100 or h < 100:          # 窗口还没完成布局（启动瞬间）
            w, h = 900, 560

        photo = getattr(self.canvas, "_tkphoto", None)
        need_rebuild = (photo is None
                        or photo.width() != w or photo.height() != h)

        if need_rebuild and hasattr(self.canvas, "resize"):
            # matplotlib 的正规入口：重建位图 + 同步 figure 尺寸
            self.canvas.resize(_SizeEvent(w, h))
        else:
            self.fig.set_size_inches(w / self.fig.dpi, h / self.fig.dpi)

    def _on_canvas_resize(self, event):
        """窗口大小变了就重画；尺寸没变（<4px）就不动，避免和重绘互相触发。"""
        if (abs(event.width - self._last_canvas_size[0]) < 4
                and abs(event.height - self._last_canvas_size[1]) < 4):
            return
        self._last_canvas_size = (event.width, event.height)
        self._redraw()

    def _redraw(self):
        orgs = self._selected_species()
        kind = self.chart_var.get()
        self._sync_figure_size()        # 先对齐尺寸（会重建位图），再清空重画
        self.fig.clear()

        try:
            if kind == "ENC-plot":
                ax = self.fig.add_subplot(111)
                plot.plot_enc(ax, self.df, orgs)
            elif kind == "RSCU 热图":
                # 热图是"宽扁"结构：单独给一条横向区域，
                # 否则 5 行会被拉得又高又怪（画布是宽屏的）
                ax = self.fig.add_axes([0.10, 0.32, 0.78, 0.48])
                sub = self.rscu.loc[[o for o in orgs if o in self.rscu.index]]
                plot.plot_rscu_heatmap(self.fig, ax, sub)
            elif kind == "PCA 散点":
                ax = self.fig.add_subplot(111)
                plot.plot_pca(ax, self.pca, orgs, show_cluster=True,
                              explained=plot.load_explained())
            else:
                ax = self.fig.add_subplot(111)
                plot.plot_box(ax, self.df, orgs)
            # 热图用的是固定位置坐标轴（add_axes），与 tight_layout 不兼容，
            # 对它调用会刷 UserWarning 且布局可能被改坏，所以跳过。
            if kind != "RSCU 热图":
                self.fig.tight_layout()
            self.canvas.draw()
            self.status.config(text=f"已绘制：{kind}　物种 {len(orgs)} 个　"
                                    f"CDS {int(self.df.organism.isin(orgs).sum())} 条")
        except Exception as exc:                      # noqa: BLE001
            self.status.config(text=f"绘图失败：{exc}")
            messagebox.showerror("绘图失败", str(exc))

    def _export(self):
        path = config.FIG_DIR / f"gui_export_{self.chart_var.get().replace(' ', '_')}.png"
        self.fig.savefig(path, dpi=220, bbox_inches="tight")
        self.status.config(text=f"已导出：{path}")
        messagebox.showinfo("导出成功", f"已保存到：\n{path}")


def main():
    app = CodonApp()
    # --selftest：自动画一遍再退出，用于检查环境（不需要人工点）
    if "--selftest" in sys.argv:
        for kind in CodonApp.CHARTS:
            app.chart_var.set(kind)
            app._redraw()
            app.update()
        print("selftest OK:", CodonApp.CHARTS)
        app.destroy()
        return
    app.mainloop()


if __name__ == "__main__":
    main()
