# 细菌密码子使用偏好分析（Codon Usage Bias）

《用 Python 玩转数据》实践项目 —— 5 种细菌全基因组 CDS 的密码子偏好量化、可视化与聚类

***

## 这个项目在做什么

**从 NCBI 取 5 个细菌的全基因组 CDS，自己算三个密码子使用指标，看不同物种的密码子偏好差多少、能不能靠密码子偏好把物种自动分开。**

四个步骤对应课程的四块内容：

| 步骤 | 脚本                      | 做什么                         | 对应课件                   |
| -- | ----------------------- | --------------------------- | ---------------------- |
| 1  | `01_fetch_data.py`      | 从 NCBI Entrez 取全基因组 CDS 并缓存 | 第 5 周 数据获取与预处理         |
| 2  | `02_compute_metrics.py` | 清洗 CDS，算 RSCU / GC3s / ENC  | 第 3–4 周 数据结构、DataFrame |
| 3  | `03_analyze.py`         | 相关性分析、PCA 降维、K-means 聚类     | 第 6 周 数据统计挖掘           |
| 4  | `04_plot.py`            | 出 4 张图（含箱线图）                | 第 4 周 箱形图、第 6 周 可视化    |
| 5  | `05_gui.py`             | tkinter 图形界面浏览结果            | 第 7 周 面向对象与 GUI        |

***

## 目录结构

```
实践项目_密码子偏好/
├── config.py                 全局配置：物种、登录号、过滤阈值、配色
├── utils.py                  核心算法：RSCU / GC3s / ENC / ENC 期望曲线
├── 01_fetch_data.py          取数（带缓存与重试）
├── 02_compute_metrics.py     清洗 + 计算指标
├── 03_analyze.py             统计 + PCA + K-means
├── 04_plot.py                出图
├── 05_gui.py                 图形界面
├── data/                     原始 CDS 缓存（*.fna，约 23 MB）
├── figures/                  输出图片（220 dpi）
├── output/                   输出表格与统计结果
└── README.md
```

***

## 环境准备

| 依赖                          | 谁需要          | 说明                           |
| --------------------------- | ------------ | ---------------------------- |
| pandas / numpy / matplotlib | 全部           | 数据分析与绘图的基础库                  |
| tkinter                     | 05（GUI）      | Python 自带，不用装                |
| scikit-learn                | 03（可选）       | 只用于 K-means；没装会自动退回 numpy 实现 |
| **biopython**               | **只有 01 需要** | 联网取数用 `Bio.Entrez`           |

也就是说：**数据已经缓存在 `data/` 的前提下，02\~05 不需要 biopython 也能跑**
（`utils.py` 内置了遗传密码表，并与 Biopython 对拍保证一致）。

```bash
# 最少装这几个就能看图
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple pandas matplotlib

# 要自己重新取数 / 用 sklearn 版 K-means，再装这两个
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple biopython scikit-learn
```

> 先改 `config.py` 里的 `NCBI_EMAIL` 为你的邮箱，NCBI 要求提供联系方式。
>
> 想确认密码表没被改坏，可以跑 `python utils.py`，它会打印"内置表与 Biopython 对拍结果"。

***

## 运行顺序

```bash
python 01_fetch_data.py --verify   # ① 先核实登录号（第一次务必跑）
python 01_fetch_data.py            # ② 下载 CDS，缓存到 data/
python 02_compute_metrics.py       # ③ 算指标
python 03_analyze.py               # ④ 统计与聚类
python 04_plot.py                  # ⑤ 出图
python 05_gui.py                   # ⑥ 打开图形界面
```

**关于联网**：只有第 ①② 步需要网。数据一旦缓存到 `data/`，后面全程离线可跑。
所以演示前先跑一次 ①②，之后即使断网也不影响。

**关于重跑**：②③④ 是纯计算，可以随便重跑，结果一样（K-means 固定了 `random_state=42`）。

***

## 输出的文件

| 文件                                 | 内容                                    |
| ---------------------------------- | ------------------------------------- |
| `output/codon_metrics_per_CDS.csv` | 每条 CDS 一行：物种、基因、长度、GC3s、ENC、59 个 RSCU |
| `output/species_summary.csv`       | 每个物种一行：CDS 数、平均 GC3s、ENC 均值/极值/标准差    |
| `output/rscu_matrix.csv`           | 物种 × 密码子的 RSCU 均值（画热图用）               |
| `output/filter_log.csv`            | 过滤统计：每种原因剔除了多少条（**报告方法节必须引用**）        |
| `output/pca_coordinates.csv`       | 每条 CDS 的 PC1/PC2 坐标 + 所属聚类            |
| `output/analysis_result.json`      | 相关系数、PCA 解释方差、聚类纯度等                   |
| `figures/fig1_ENC_plot.png`        | GC3s–ENC 散点 + 无选择期望曲线（核心图）            |
| `figures/fig2_RSCU_heatmap.png`    | RSCU 热图（按氨基酸分组）                       |
| `figures/fig3_PCA.png`             | PCA 散点（按物种着色）                         |
| `figures/fig4_GC3s_box.png`        | GC3s 箱线图                              |

***

## 常见问题

**Q1. 报 `Failed to understand id: GCF_...`**
NCBI 的 nuccore 库不认 GCF\_ 组装号，必须用染色体自身的 `NC_` 登录号。本项目 `config.py` 里用的都是 NC\_ 号，如果你自己加物种，记得换。

**Q2. 报连接被重置 / 超时**
NCBI 偶发限流。脚本已内置重试（`config.NCBI_RETRY`）。连续失败就过几分钟再跑，或把 `config.NCBI_SLEEP` 调大一点。

**Q3. 图上中文显示成方框**
`utils.setup_chinese_font()` 会自动找 SimHei / 微软雅黑。Linux 上要自己装中文字体。

**Q4. 没装 scikit-learn 能跑吗**
能。`03_analyze.py` 会自动退回自带的 numpy 版 K-means，只是会慢一点。PCA 本来就是手写 SVD，不依赖 sklearn。

**Q5. 想换物种**
改 `config.py` 的 `SPECIES` 字典（加 `NC_` 登录号 + 遗传密码表编号），再跑一次 `01_fetch_data.py --verify` 核对，然后按顺序重跑 ①\~④。

> ⚠️ 遗传密码表别乱填：细菌一般用 `11`，但支原体要用 `4`，线粒体用 `2`。
> 填错了所有指标都会算错，而且不会报错。

**Q6. 报 `UnicodeEncodeError: 'gbk' codec can't encode character ...`**
中文 Windows 的控制台是 GBK 编码，**`print()` 里不能出现 emoji**（⚠️、✅ 之类），否则脚本会在打印那一行时直接崩掉（哪怕前面的计算都成功了）。本项目所有 `print` 已改用 `[注意]`、`[完成]` 这类纯中文标记。自己加提示语时注意这一点。

**Q7. 看到 `can't open file for write !`**
这是 matplotlib 写字体缓存（`~/.matplotlib`）时的提示，不影响出图。若图确实没生成，检查 `figures/` 目录是否存在或是否只读。

***

## 复现说明

- 随机种子：`04_plot.py` 抽样用 `seed=42`，`03_analyze.py` K-means 用 `random_state=42`
- 数据版本：以 `data/` 下缓存的 FASTA 为准（NCBI 记录会更新，所以缓存必须保留）
- 过滤参数集中在 `config.py`，报告里要写清用的是哪一套参数
