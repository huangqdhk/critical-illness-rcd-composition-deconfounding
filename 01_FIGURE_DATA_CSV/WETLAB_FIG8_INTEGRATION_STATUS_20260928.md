# 湿实验 Figure 8 接入图数据包 · 现状与待办（2026-09-28）

> 本文件是**核查台账**，不是面板数据。文件名 `WETLAB_FIG8_INTEGRATION_STATUS_*` 不匹配 `Figure_*` 面板命名，不会被 `FIGURE_CROSSWALK.csv` / manifest / lint 当作面板文件。
> **本次核查未改动、未改名、未移动 `01_FIGURE_DATA_CSV` 下任何文件。**

---

## 一、一句话结论

湿实验 Figure 8 **尚未接入** `01_FIGURE_DATA_CSV`（交叉表里只有 Figure 1–7 + S1–S14）。接入前有**一个必须先解决的静默覆盖陷阱**（第二节），且最终面板集取决于导师删 W2/W3 的结果，**两项都需等 v6 与新 Figure 8 成图**。

---

## 二、⚠️ 静默覆盖陷阱（必须先解决，否则会污染已发布源数据）

`Figure_8A.csv` … `Figure_8G.csv` 这 **7 个文件名在两个包里同时存在，但内容完全不同**：

| 文件名 | `01_FIGURE_DATA_CSV/Main/`（pkg 世代 = **正文 Figure 6**） | `E:\SCI\Figure_8\Figure_8\`（= **湿实验 Figure 8**） | 同名同内容？ |
|---|---|---|---|
| Figure_8A.csv | 812 B — 表头 `cohort,n_obs,n_patients,beta1_std,se,p_one_neg,...` | 13,102 B — 表头 `实验,样本,组,基因,臂,Cq均值,内参Cq(Gapdh),...` | **不同** |
| Figure_8B.csv | 529 B | 9,678 B | **不同** |
| Figure_8C.csv | 937 B | 6,796 B | **不同** |
| Figure_8D.csv | 729 B | 8,512 B | **不同** |
| Figure_8E.csv | 386 B | 5,056 B | **不同** |
| Figure_8F.csv | 1,279 B | 3,344 B | **不同** |
| Figure_8G.csv | 3,154 B | 34,407 B | **不同** |
| Figure_8H.csv | *不存在* | 17,114 B | — |

**含义**：包内 `Main/Figure_8A–8G.csv` 是**正文 Figure 6**（时序与干预响应；manifest 指向 `_intermediate/M11M12_step*`，图注 8E/8F/8G 分别对应正文 6E/6F/6G）。湿实验 Figure 8 的源包 `E:\SCI\Figure_8\Figure_8\` 用了**同名文件**装**完全不同的数据**。

> **因此：绝不能把 `E:\SCI\Figure_8\Figure_8\*.csv` 直接复制进 `01_FIGURE_DATA_CSV/Main/`** —— 会静默覆盖正文 Figure 6 的 7 个面板源数据，并顺着 `FIGURE_CROSSWALK.csv` → `N5b_build_source_data_xlsx_crosswalk.py` 污染 Source Data 的 xlsx 与投稿包。这正是本包 README 反复警示的"手工副本必然漂移"类事故。

### 需要定的命名口径（三选一，建议 ①）

1. **给湿实验 Figure 8 一套独立 pkg 名**（如 `Figure_8w A–H` 或 `Figure_10 A–H`），在 `FIGURE_CROSSWALK.csv` 里以 `ms_figure=Figure_8` 显式桥接。与包内"pkg 世代 ≠ 正文图号"的既有约定一致，不动任何旧文件。
2. 把包内旧 `Figure_8A–8G.csv` 改名为其对应的正文图号（`Figure_6A–6G`），腾出 `Figure_8*` 给湿实验图。—— 但会**打断 pkg 世代命名惯例**，且 `FIGURE_CROSSWALK.csv`、`FIGURE_DATA_MANIFEST.csv`、Source Data 生成器全需同步改，风险最大。
3. 湿实验 Figure 8 数据**不入本包**，仍以 `E:\SCI\Figure_8` 为唯一源，只在 `FIGURE_CROSSWALK.csv` 里登记 `csv_relpath` 指向外部路径。—— 最省事，但破坏"包自包含"约定，Source Data xlsx 生成器需支持外部路径。

---

## 三、其它两项待决

**1. 湿实验 Figure 8 的最终面板集未定**。盘上新 Figure 8 成图（`可视化\01_主图\Figure_8`，2026-09-28 21:50 重建）与图注（`可视化\03_审计台账\Figure_8.legend.md`）目前均为 **8 面板 a–h**：a–e 用 W1、**f = W2**、**g = W3**、h = THP-1。导师"删除了 W2/W3 的相关内容"若落到 Figure 8 上即撤 f、g 两面板（剩 a b c d e h 六个），但**是否撤、撤后是否重排字母（h→f）尚未确认**。字母重排会连带改全部面板文件名与图注。

**2. 陈旧遗留 `Main/Figure_9A–F.csv`（6 个文件，共 22,991 B）**。按 `CDD投稿需要.md` 记录，它们是"旧 9 图方案遗留（内容 = 重编号前的旧 Figure 10/11）"，终版无 Figure 9，需确认去向（转 SI 数据源或移交 `归档/`），**勿混入投稿包**。本次未改动，仅登记。

---

## 四、v6 与新 Figure 8 到齐后的执行清单

1. 确认最终面板集与字母顺序；核对 `Figure_8.legend.md`、成图、面板 CSV 三者一致。
2. 按第二节定下的命名口径，在本包建立湿实验 Figure 8 的面板数据（**新增文件，不改旧文件**）。
3. 更新 `FIGURE_CROSSWALK.csv`：新增 `ms_figure=Figure_8` 的行；主图由 7 张 → 8 张。
4. 更新 `FIGURE_DATA_MANIFEST.csv`：新增 Figure_8 各行（面板 → 源文件）。
5. 更新 `README_图数据包说明.md`：主图 7 张（54 面板）→ 8 张；总面板 149 → 新值；世代沿革补一条（2026-09-28 接入湿实验 Figure 8）。
6. 处置 `Main/Figure_9A–F.csv`（移入 `归档/` 并登记，或确认转为 SI 数据源）。
7. 跑 `python lint_package.py` 要求全绿；重建 Source Data xlsx（`N5b_build_source_data_xlsx_crosswalk.py`）并跑 `F0_FIGURES/sync_submission_data_copy.py --apply` 同步投稿树镜像（该脚本会以 sha256 断言两树逐字节一致）。
8. 与 `可视化` 侧成图一致性核对：面板数、N、误差线/无误差线声明（CDD 对 N 与误差棒描述均为 must）。

---

## 五、本次未做的事（明确声明）

- 未向本包新增/改名/移动/删除任何文件。
- 未改动 `FIGURE_CROSSWALK.csv`、`FIGURE_DATA_MANIFEST.csv`、`README_图数据包说明.md` 及任何 `Main/`、`Supplementary/` 面板 CSV。
- 未触碰 `04_AUDIT_GOVERNANCE/`、未改动主文稿。
