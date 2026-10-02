# RESULTS_MANIFEST_v1.0 — 结果表包冻结说明（D4/N22 治理）

**2026-09-21 归档区分类整理（登记层变更，不动任何数据/报告）**：`归档/` 由「13 个平铺子目录 + 16 个顶层散件」重整为 **9 个类别目录**（`01_第1轮R1批次` / `02_模块调试批次` / `03_一次性脚本与运维` / `04_文稿与计划备份` / `05_清单与图数据备份` / `06_工具脚本备份` / `07_审计报告备份` / `08_投稿树旧世代与备份_20260921` / `09_期刊指南文本提取_202609`），细目见 `归档/README_归档说明.md`。本 manifest 同步：**306 行 archived 的 `location` 列改写**（MISSING_OK 两行 location 保持空、`README_归档说明.md` 保持 `归档`）+ **新登记 8 件**（`_cdd_gta_text.txt`、`_aj_artwork_text.txt`、`_dump_gta.py`（自 `E:\SCI\` 根入档）、`RESULTS_MANIFEST_v2.0_old11fig_backup.csv`、`Cover_Letter_CDD_draft_20260905.md`（自包根入档）、`README_期刊指南提取说明.md`（新建）、以及本批次两份执行留痕脚本 `归档分类整理_20260921.py` / `归档分类整理_manifest同步_20260921.py`）+ 本行随 `README_归档说明.md` 重写刷新 size/sha。归档区与登记仍按 basename 一一对应（315 = 315），无重复 filename、无 canonical/archived 同名冲突。因归档件位移而同步路径的现役引用：`P0_manuscript_source_check.py`、`P0_manuscript_source_check_ext.py`、`P13_crossref_figure_check.py`（读归档中文稿）与 `03_LOGS/_diag_{final,md,rowcompare}.py`、`_gen_teacher_csv.py`（读旧 11 图快照）。写入备份：`_review_tmp/RESULTS_MANIFEST_v2.0.bak_before_归档分类_20260921.csv`。

**2026-09-21 figure_data 段重登记（原审稿意见 §9-2 第三项收口；原文见 归档\07_审计报告备份\可视化文档_合并前_20260928/）**：图源数据树 `01_FIGURE_DATA_CSV/{Main,Supplementary}` 的文件名沿用 **pkg 世代**编号（09-05 九图体系，如 Figure_3C/5J/8A/S9F），正文引用用 **ms 世代**编号（Figure 2G/4J/6A/S11F），两者由 `可视化/00_对照表/FIGURE_CROSSWALK.csv`（148 行，ms_panel ↔ pkg_panel ↔ csv_relpath）桥接。本段此前只做了一半迁移：19 行改成了 ms 名（磁盘上不存在 → 幽灵行），其余仍是 pkg 名但 size/sha/行数多为迁移前旧值（96 行"同名不同面板"），于是 lint 的 SHA256 校验长期在拿错位的两份文件互比，**29 个在用图源文件根本没进校验圈**。本次以 crosswalk × 磁盘真值重建登记：**149 个图源 CSV 与 149 条 canonical 行一一对应**（刷新 120 / 改名 10 / 新增 19），19 条 ms 名幽灵行转 `status=superseded` 仅作历史留痕（lint 不再校验，location 置空以免 fp() 指向错位路径）。**重登记只改登记、不动数据**：无任何 CSV、图件或数值被修改。重登记前备份：`归档/05_清单与图数据备份/RESULTS_MANIFEST_v2.0_pre_figure_data重登记_20260921.csv`（sha256 见其登记行；影子备份在 `_review_tmp/RESULTS_MANIFEST_v2.0.pre_before_fix92.csv`）。写入纪律：先备份 → 写 `.tmp` → 复核行数与"canonical 图源行 == 磁盘 CSV"不变式 → `os.replace` 原子替换（本批次第一次写入曾因下述坏行触发 `writerows` 异常，因该纪律未造成损失）。

**2026-09-21 同日修缮四项**：① **P1_IIAMD_signature_v1.0.csv 登记行残余缺陷补修**：09-16 记录所称的该行修复只落到数值列（size_bytes/n_data_rows/sha256 已就位），`analysis` 字段**仍缺引号**，字段内逗号把整行切成 16 列 → `csv.DictReader` 会截断 analysis、把 notes 错置为 "367/down 4"，多余值落到 `None` 键；本次按字段语义回接（前 11 列原样 + analysis 三片用逗号回接 + notes + 空 class_），行首次可正确解析为 14 列，逐字符保留原文。② **lint 新增 `read_data_header()`**：跳过**全部**前导 `#` 注释行（原实现只跳 1 行；m5 拆出的 `Figure_S2E_fit_stats.csv` 带 2 行溯源注释，被误判"缺 gene_set_version/score_version 列"，其实两列俱在且值为 NA）。③ **FDR 空值豁免键随重登记改为磁盘真名**：`Figure_5F.csv`/`Figure_5G.csv`/`Figure_7E.csv`/`Figure_S6I.csv`（正文名 Figure 4F/4G/S10D/S7I 写进理由文本）——原键名是 09-16 那轮改成的正文图号，与磁盘名错位后豁免失效。④ **版本列豁免清单体检**：移出 13 个"磁盘文件本身带版本列"的条目（Figure_5E–5I / Figure_6A–6F / Figure_S10C / Figure_S10D，豁免等于把可校验文件排除在校验外）与 1 个死条目（Figure_6G.csv，磁盘无此文件）；补入 11 个此前从未被检查过的模块图数据文件（正文 Figure 6A–6G、7A–7C、S10C 对应的磁盘名 Figure_8A–8G、9A–9C、7G），按既有 M10–M16"导出设计不含版本列"政策逐条取证登记。**收口实测：lint 由 通过 73 / 失败 162 → 通过 76 / 警告 0 / 失败 0；受 SHA256+尺寸校验的 canonical 文件 372 → 500；受版本戳校验的表 272 → 359。** 残留（自披露，不属本批次）：11 个模块图数据文件本身无版本列，要补列须重跑对应导出脚本（产物改动，本次未做）。

**2026-09-16 A+B 主图压缩批次（迎接湿实验主图 Figure 8）**：主图 9→7（54 面板：6-8-8-10-9-7-6）、附图 12→14（94 面板）。① 合并：旧 Figure 2+3 → 新 Figure 2（模块级解离＋跨平台/外部验证，8 面板：旧 2A–2D＋旧 3A/3B/3C/3E）；旧 Figure 6+7 → 新 Figure 5（统计单位重算＋组成审计＋扰动三角，9 面板：旧 6A/6D/6E/6F＋旧 7A–7D/7F）；旧 Figure 4→3、5→4、8→6、9→7 顺移。② 降级：10 个面板降入附图——旧 2E/2F/3D/3F/3G/3H → 新 Figure S2A–S2F（§2/§3 首引）；旧 6B/6C/7G/7E → 新 Figure S10A–S10D（§10/§11 首引）。③ 附图级联（保持全文引用严格升序，同 2026-09-05 二次重排口径）：S2→S3、S3→S4、S4→S5、S5→S6、S6→S7、S7→S8、S8→S9、S9→S11、S10→S12、S11→S13、S12→S14。④ 面板 CSV 两阶段安全重命名 133 个（主图 54＋附图 79），**数值零改动、sha256 不变**，映射台账 `01_FIGURE_DATA_CSV/FIGURE_RENUMBER_MAP_20260916.csv`；RESULTS_MANIFEST 133 行行级 filename/location 同步（非 pandas 往返，规避 B1 批次 NA 误伤模式），磁盘↔manifest 双向核验 148/148 零缺失、sha 抽查随名迁移一致。⑤ 文稿 `111文稿_v5_英文版.md` 全文图引用重映射（67 项复合引用占位替换＋单令牌正则重映射，裸图号 0 处）＋图例块重写（新 F2/F5 合并图注、新 S2/S10 图注、组图注行与 causal note 顺移）；**顺带修复一处 v5 遗留误引**：Discussion "(Section 7, Figure 5B)"（尸检肺铁轴蛋白）实指铁轴蛋白面板 → 改为新 Figure 4G（原应引 5G）。⑥ lint 同步（M1X 白名单 9D-F→7D-F、7A-G→5E-I/S10C/S10D、8A-G→6A-G、9A-C→7A-C；FDR_NULL_OK 键 S6I→S7I、5F→4F、5G→4G、7E→S10D；M2 守门 Figure_5A-E→4A-E）；P13 交叉引用核查脚本期望值同步。⑦ **遗留待 Python 环境执行**（本机无 Python/R，所有 .py 改动均未运行）：`lint_package.py` 复跑须 74/0 全绿；`投稿CELL DEATH AND DIFFERENTIATION/Source_Data/` 22 个 xlsx 须由 `N5_build_source_data_xlsx_20260905.py` 重新生成（现仍为旧 9 图体系）；`P13_crossref_figure_check.py` 复跑出报告。⑧ 不同步项（自披露）：归档/04_文稿与计划备份/111文稿_v5.md 中文版与"-手动步骤的修复"兄弟夹仍为旧图号（归档只读）；05_RELEASE_GITHUB/repo 为 2026-08 快照旧编号（定稿发布时重切）；湿实验 Figure 8 与新 S15/S16/S17（含 WB 全膜图）面板位已预留，待真实湿实验数据审计通过后并入（见 `111湿实验/08_WB到达后_主图附图设计方案_20260916.md`）。改动前备份：归档/ 下 `*_pre_AB压缩备份_20260916.*`（文稿/manifest/lint/图数据 zip，2026-09-21 分类后位置见 `归档/README_归档说明.md`）。

**同日复跑修补记录（2026-09-16 晚，Python 3.14.5 @ D:\MovedFromC\Python314）**：lint 复跑暴露并处置四项——① **P1_IIAMD_signature_v1.0.csv 登记行历史缺陷修复**：score_version 列含未加引号逗号（"up 4,367/down 4,720"）致整行右移、size_bytes 错挂 "720）"，致 lint 在 checksum 循环即崩溃（该缺陷在 09-16 上午备份中同样存在，属 v2.0 重建遗留）；已按磁盘实况（size 2,474,865、sha256 不变）行级修正列归位，score_version 记文件真实戳。② **S93/S93b 白名单补登**：lint 的 M1X_NO_VERSION_STAMP 缺 M17 批次两表条目（09-05 变更记录称已入但实际缺失），已补登。③ **百度云同步残留清理**：04_AUDIT_GOVERNANCE 下 7 个 `.baiduyun.uploading.cfg` 隐藏残留文件删除。④ **归档备份登记补齐**：归档/ 下 11 个未登记备份文件（本批次 4 个 + 09-05/06/07 会话遗留 7 个）行级登记为 archived 行（sha256/size 按本机磁盘实况）。**遗留环境项（自披露，不属本批次）**：本机 归档/ 为百度云部分同步状态，登记在册的 ~300 个历史归档文件本地不存在，致 lint 两项"登记在册但本地不存在/归档不一致"在本机无法转绿——须在同步完成的主机复跑确认；本机复跑结果 72 通过 / 0 警告 / 3 失败（失败均为该环境项与本文件自身校验和滞后，后者随本行登记即时更新）。

**2026-09-06 遗留问题修复批次（N5c）**：对 2026-09-05 质控报告的四项作者侧遗留逐项处置——① **中文版 v5 图号同步**（归档/04_文稿与计划备份/111文稿_v5.md，改动前备份 111文稿_v5_pre_图号同步备份_20260906.md）：双遍令牌重映射（图/Figure 双前缀，121 处引用列表段零丢失）+ 图例区整段重写（镜像英文终版结构，含 Figure 5B 图注 k=3/k=16 纠偏）+ 8 处引用补充移植（Figure 1B/1C 与 S1–S12 全集首引）+ 顶部同步声明 banner；正文被引 89 面板全部对得上磁盘；该文件不在冻结链（无 sha 登记），banner 内自披露。② **兄弟夹"-手动步骤的修复"处置**：新增 文件夹状态声明_20260906.md（历史副本+旧图号警示+恢复基底不可删）；修复其 CRediT 草案 "Figure 1–11"→"Figure 1–9"（附重排说明）。③ **Supplementary Notes 披露**：新增 README_SupplementaryNotes_20260906.md（NicheNet 报告 Figure_13* 远古编号遗物披露，不改写冻结文件）。④ **Table S93/S93b 内容审计**（M17_S93_Data_Audit_20260906.md）：13/13 检出基因统计全复现（Welch 口径=log2 变换后，Δ≤5×10⁻⁴；首轮原始尺度复算偏差系审计口径误设，已锁定披露）、4 缺行基因=未检出设计内缺席、PXD050432 外部活验=Front Pharmacol 2024;16:1397498 异去氢钩藤碱(LPS_IRN) 小鼠肺 ALI 蛋白组——**审计通过**；遗留：S93 尚未被文稿引用（M17 会话入稿方式待定）+ 建议表注补 log2-Welch 口径句。**另发现并修复 A 类一处**：两英文稿作者行仍为 Huang–Wang–Li，与已落稿等贡献声明（Huang–Li–Wang）及 CRediT 状态注记矛盾——按导师确认口径补完两夹作者行同步（仅换位，各自上标不变）。新增 2 个治理文件入册（M17 审计报告 + Notes README）；**lint 全绿：74 项通过 / 0 失败**（2026-09-06）。

**2026-09-05（下午）严格质控 + 附图二次重排批次（N5b）**：对上午 CDD 重排交付执行主编级全面质控，发现并修复五项：① **附图首引顺序违规**——上午补入 S1–S7 正文首引后附图首引不再升序（Nature Portfolio 要求附图按引用顺序编号），按首引位置二次重编号（79 个 CSV：S7→S2、S2→S3、S8→S4、S4→S5、S9→S6、S10→S7、S5→S8、S6→S9、S11→S10、S12→S11、S3→S12），manifest 79 行同步、lint 白名单键 S9I→S6I、文稿附图令牌整体重映射 + 散文段/causal note 人工重写，主图编号不动；② **图注-数据错配（v5 原有缺陷）**——正文与图注将 k=16 池化值（+1.50/−0.94/+1.02，实源于 Figure_6A 层）挂于 Figure 5B 名下，而 5B 面板数据实为六队列逐层效应+历史 k=3 池化行（+1.85）；已改 §6 引用指向 Figure 6A 并重写 5B 图注如实描述；③ **Figure 1B/1C 与附图 S1–S7 正文从未被引（v5 原有）**——已补引用；④ **Excel 源数据两处保真缺陷**——布尔列（heidi_pass/same_sign/sex_included）Excel 往返漂移为 1.0/0.0、Figure_S12D（旧编号）注释首行被误当表文致 sheet 损坏；生成器已修（布尔文本保真+注释剥离入 Panel_X_notes），全量 148 面板复验零差异；⑤ **脚本再生旧名风险**——历史导出脚本仍写旧图号文件，按任务8 先例保留原样，README 已加显式警示。**事故披露**：质控中 FIGURE_DATA_MANIFEST.csv 曾因写回异常被截断，已以兄弟交付夹（-手动步骤的修复/01_FIGURE_DATA_CSV）144 行副本为基底、经两次映射+6 面板补登完整重建（149 行，含双迁移链注记）；该文件未入 RESULTS_MANIFEST 冻结链（无 sha 登记），此前亦无既定备份纪律，本次起随包维护。数值层零改动（全部 CSV sha256 不变）；**lint 全绿：74 项通过 / 0 失败**（2026-09-05 下午，含 S93 补登记）。遗留作者侧：中文版 v5（归档/）与兄弟夹"-手动步骤的修复"内英文版副本仍为旧图号（未同步，待作者裁定）；Supplementary Notes 中 NicheNet 报告含远古 Figure_13* 引用（历史文件，不影响投稿包）。

**2026-09-05 CDD 投稿适配图号重排批次（N5）**：为匹配 Cell Death & Differentiation 惯例（主图每图 ≤8–10 面板、各图面板数均衡、全文引用严格升序），主图 11→9（64 面板：6-6-8-8-10-6-7-7-6）、附图 12→12（84 面板：5-10-8-9-6-8-7-5-9-7-5-5）。① 面板 CSV 两阶段安全重命名 91 个（主图 41 + 附图 50；映射台账 `01_FIGURE_DATA_CSV/FIGURE_RENUMBER_MAP_20260905.csv`，脚本 `N5_figure_renumber_20260905.py`）——**数值零改动，sha256 不变**，manifest 91 行 filename 同步并附迁移注记；② 合并逻辑：旧图 4+5→新图 5（Part II 临床+生化边界）、旧图 8+9→新图 8（时序+干预）、旧图 10+11→新图 9（工具+力学）；附图旧 S3(EFA) 并入新 S2、旧 S12 拆分（因果层 5 面板→新 S8、ML/DL/表观 5 面板→新 S12）；③ 消除旧图 2 的 16 面板超载、旧图 9/10/11 各 3 面板过薄、"§2 引图 3、§3 引图 2I–2P"乱序三缺陷，并补齐旧图 6 正文引用的历史缺口；④ lint 同步（白名单/FDR_NULL_OK/M2 守门内图号迁移注记），**lint 全绿：74 项通过 / 0 失败**（2026-09-05）；⑤ `FIGURE_DATA_MANIFEST.csv` 重建（149 行，补登 M15/M16 六面板 9A–9F），README 图数据包说明同步；⑥ 投稿图源数据按 Nature Portfolio 惯例每图合并为单 xlsx（每面板一 sheet，22 个工作簿，`投稿CELL DEATH AND DIFFERENTIATION/Source_Data/`，脚本 `N5_build_source_data_xlsx_20260905.py`，抽查数值零差异）；⑦ 登记补录：当日上午另一会话 M17 批次产物 Table_S93/S93b（PXD050432 急性 LPS 肺蛋白组 WB 靶标核验）登记入册并入 lint 白名单。文稿 111文稿_v5_英文版.md 全文图引用与图例块同步更新（§2/§3/§4/§6/§7/§8/§10/§13/§14/§15/§16 + 图例 1–9/S1–S12）。

**2026-08-29 DOI 回填与注册 ID 纠错（投稿准备批次）**：三份 OSF 注册已于 2026-08-29 解禁公开并取得 DOI——**C7RYD**（冻结分析计划，10.17605/OSF.IO/C7RYD，2026-08-20 注册）、**ETVMJ**（M10/M13/M14 二次数据预注册，10.17605/OSF.IO/ETVMJ，2026-08-24 注册）、**98CM3**（M11/M12 预注册，10.17605/OSF.IO/98CM3，2026-08-27 13:45 注册，OSF 官方 date_registered=2026-08-27T05:45:15Z=北京时间 13:45:15）。DOI 已回填四处：注册文档头（M10_M13_M14_pre_registration_20260824.md / M11_M12_pre_registration_20260827_draft.md）、双版文稿注册表述（版本说明、4.22、4.25）、本文件变更行、冻结计划 DOI 占位行（P0_FROZEN_ANALYSIS_PLAN_v1.0.md）。本地文档误记的注册 ID（3y2rb、3e92j）已全局纠正为真实 ID（ETVMJ、C7RYD；39 文件 108 处，替换日志 _intermediate/DOI_replace_log.txt）；OSF API 核验 date_registered=2026-08-27T05:45:15Z = 北京时间 13:45:15，A1 时序裁定获官方记录二次确认。2026-08-29 终稿按作者指示：双版文稿内注册号（osf.io/xxx）全部移除、仅保留 DOI 表述。

**2026-08-29 B1 溯源修订（v4.7 终审未竟项批次）**：`M13_Confounding_Audit_Report.md` 数值修订（§1 补录 myeloid 零模型 bootstrap 数值对：null beta_myeloid=1.612±0.063、经验 p=0.0005；来源 `M13_step2_null_bootstrap.py` B=2000 seed=0，B=10,000 稳健核验 p=0.0001/零分布 1.613±0.062），manifest 行 SHA256/size/n_data_rows 同步更新；新增 `03_LOGS/M13_null_bootstrap_myeloid_log.txt`（核验日志，03_LOGS 不纳入 manifest 登记）。**事故记录**：本次修订过程中 manifest CSV 曾因 pandas 往返写回将字面量 "NA" 单元格误写为空串（并引入浮点化/换行符差异），已按 `P0_rebuild_manifest_v2.py`（v1.0 基底 + 磁盘实况）重建参考文件逐格比对修复（恢复 955 格），notes/analysis 治理文本保留；同时恢复了 M1x_step12 清理批次 64 行的登记约定（3 个更新行 location/status=archived + 61 个新增行 class/版本列空串），lint 复绿 74/0。

**2026-08-21 v2.1 变更记录（投稿前硬阻断项收尾批次）**：RESULTS_MANIFEST_v2.0.csv 以 v1.0（279 行，2026-08-17 冻结）为基底重建为 591 行（重建脚本 `P0_rebuild_manifest_v2.py`；v2.0 曾因修复脚本写回异常被截断，重建以磁盘为权威重算全部 canonical 文件 SHA256，v1.0 文件保留历史冻结链）。本批次变更：① 修正 M2/M3/M4 批次引入的版本戳漂移（manifest 空串 vs 冻结 CSV 'NA'，164 行）；② 5 个 CSV 补版本列（Figure_10D、S15e Instruments/LeaveOneOut/Sensitivity、S60）+ 4 个 Figure CSV（1B/1C/2E/S2F）；③ N22 s_number 归一（S15e/f/g→S15、S19i→S19）；④ archive 子目录重组后 location 修正（R1_table_backups/M1_check_debug/R1_artifacts）；⑤ 补登记 P0/P1/P2 系列与 2026-08-21 新增治理文件（P0-8、查新 v1.1/重跑日志、参考文献核实报告、manuscript-to-source 报告、软件快照、pip freeze×2、P1-4 报告等）；⑥ lint 同步更新（FDR_NULL_OK 增 M4 蛋白层 4 表与 Figure_5I；版本戳多值 ';' 集合语义；archive 递归比对）；**lint 全绿：74 项通过 / 0 失败**（2026-08-21）。

**2026-08-31 M16 力学边界模块批次**：新增第九条候选边界（力学边界）全套产物——5 个附表（S88 mech_v1.0 模块 23 基因三层+PMID 活验溯源 / S89 四套 VILI 同口径对比+合并 / S90 GSE2411 力学×LPS 2×2 互作 / S91 人体桥接+Visium 空间 / S92 同源映射覆盖审计）+ 3 个主图数据（Figure_11A/B/C）+ 2 个治理文件（GEO 核验、执行报告）；预注册 M16_pre_registration_20260831.md 按 M1/M2/M10 纪律不入册（根目录治理文档）；8 个 CSV 按 M10–M15 导出设计登记为无内嵌版本列（manifest 记 NA，lint 白名单 M1X_NO_VERSION_STAMP 同步补登）；双版文稿（v5 中文/英文版）同步更新：数据集表 +4 行、方法 4.28、结果第 17 节、讨论四向→五向边界、图 11 注、参考文献 59–66（Vancouver，PMID 全经 eutils 活验）、数据可得性。**lint 全绿：74 项通过 / 0 失败**（2026-08-31）。判据裁定：A 阴性、B 形式未达（直接增强显著）、C 阴性/混合、D 达提示级——阴性/部分结果如实界定，不升级措辞。

**变更记录（References 落盘批次，2026-08-29）**：本批次变更：① P0_References_1to64_Verification_Report.md 修订（11 条待确认条目裁定记录 + round5 系列检索留痕 + GB/T→Vancouver 落盘与 README 退役记录），manifest 同步 size/sha256（8861→16653）；② 补登记 2 个治理文件：P0_Manuscript_Source_Check_Report_v2.md（§四-2 扩展终检 13/13 PASS，含 References 完整性断言）、P0_Novelty_Search_Rerun_Log_20260829.md（§四-1 查新重跑，空白维持）；修订方式为行级精确写入（非 pandas 往返，规避 B1 批次误伤模式）。**lint 全绿：74 项通过 / 0 失败**（2026-08-29）。

**冻结日期**：2026-08-15
**机器可读清单**：`RESULTS_MANIFEST_v1.0.csv`（222 行 = 211 canonical + 11 archived；含 SHA256 校验和；2026-08-15 N3 修复后重冻结新增 S2b/S8b/N3 审计文件，同日补恢复 R1 期间丢失的 `Table_S5_PCD_Scores_RAW.csv`（82 行 ARDS→Sepsis_COVID 与主表同步，数值零改动）；**同日 N2 臂定义统一后再次重冻结**——S5b 增 `arm` 列、S5 方法注双臂归组更正，2 个 `*_pre_N2_backup.*` 入 archive/；**同日 R7 结构注释更正后再冻结为 223 行 = 211 canonical + 12 archived**（`Table_S36c` PDB 注释按 RCSB 官方订正，备份入 archive/，见 §七c）；**2026-08-16 P2 批次（N3 残留重算）后冻结 233 行 = 211 canonical + 22 archived**；**同日 P3 批次（N3 残留4）后再冻结 235 行 = 213 canonical + 22 archived**——新增 `SAMPLE_MANIFEST_v1.0.csv`（样本级 manifest，800 行 = 780 sample + 7 dataset_summary + 13 external_source，见 `N3_Sample_Manifest_Report.md`）与该报告；lint 同批新增样本级不变量检查）；**同日 D1b 批次（GSE32707 外部验证，S50–S54）后冻结 242 行 = 220 canonical + 22 archived**（同载 `03_LOGS/D4_build_log.txt`）；**同日 S6 旧编号 PNG 归档后再冻结 242 行 = 219 canonical + 23 archived**——`Figure_S6_Dissociation_hexbin.png`（旧 S6 编号遗物，解离检验图现编号 S48，见 `01_FIGURE_DATA_CSV/README_图数据包说明.md`）移入 archive/，数值零改动）；**同日目录重组（导师版结构）后再冻结**——`01_RESULTS_TABLES/` 拆分为 `02_SUPPLEMENTARY_TABLES/`（`SUPPLEMENTARY_Tables_CSV` + `Supplementary_Notes`）与 `04_AUDIT_GOVERNANCE/`，`archive/` 上移根目录，日志目录 `04_LOGS`→`03_LOGS`；manifest 新增 `location` 列；`Table_S5_PCD_Scores_RAW.csv` 原文件未随包迁移、由 canonical 脚本 `build_Table_S5.py` 重生成并按 §八 重冻结（z 值逐样本一致 392/392，见该行 notes）；`Table_S18` 三表自 archive 字节级恢复（SHA 核对一致）；`Table_S1` prev38 ×2 与 `Table_S11` 大矩阵 ×2 登记缺席留痕（与 lint `MISSING_OK` 白名单同源）
**校验入口**：项目根目录 `python lint_package.py` → `03_LOGS/lint_report.md`（当前 **全绿**：60 项通过 / 0 失败）
**构建脚本**：`build_D4_package.py`（幂等，可重跑；操作日志 `03_LOGS/D4_build_log.txt`）

对应质控报告条目：**D4 治理工程化**（canonical 基因 manifest + 全表版本列 + lint）与 **N22 补充表编号冲突与版本并存**。

---

## 一、N22 编号去重决策（每个 S 编号唯一对应一个分析族）

| 冲突 | 决策 | 理由 |
|---|---|---|
| S6 双占用（Scissor vs 单细胞解离检验） | **Scissor 保持 S6**；解离检验改 **S48**（`Table_S48_Dissociation_singlecell_correlation.csv` + `Table_S48_Dissociation_Method_Note.txt`） | 方案 §8.1 原始登记 S6=Scissor（步骤 6、主图 Figure 7）；解离检验为 2026-08-13 新增分析，顺延取新号。文稿 2 处引用同步改 S48（含 Figure S6→S48） |
| S3 双占用（WGCNA 模块 vs 细胞类型组成） | **细胞类型组成保持 S3**；WGCNA 模块改 **S49**（`Table_S49_WGCNA_Modules.csv`） | 文稿 line 139 的 "Table S3" 指细胞组成（引用稳定性优先）；WGCNA 表未被以编号引用，零引用风险 |
| `Table_s12_*`（小写 s）与 S12 并存 | 改 **`Table_S12b_Mitoxyperilysis_Upstream_TF.csv`** | 并入 S12 子表序列（main + Summary + b），消除大小写双命名空间 |
| S26b 双占用（Cluster_Statistics vs HLCA 注释） | **HLCA 注释改 S26a**（方案 §26 行内既定建议）；Cluster_Statistics 保持 S26b | S26a 空闲；b–k 序列不变 |
| S18 `_CORRECTED` 与原版并存 | **`_CORRECTED` 晋升为正式表名**（Main/S18b/S18c）；原非修正版（内容已逐字节一致，2026-08-15 修正时同步）移入 `archive/*_superseded.csv` | 方案 §8.1 既定"以 CORRECTED 为准"；消除投稿包内版本并存 |
| 备份文件与顶层混杂 | 6 个 `*_pre_R1/R2/N1/fix_backup*` 全部移入 `archive/` | 备份不进入投稿包命名空间；审计链保留 |

**编号总况**：S1–S9、S11–S24、S26–S33、S36–S41、S44–S49 共 43 个在用编号 + S10/S25（跳过：空间转录组/影像组学无数据）+ S34/S35/S42/S43（未启用）。多对一映射（合法）：`Table_S22_S23_*`、`Table_S29_S31_*`（报告）、`Table_S32_S33_*`（方法注）、`Table_S38_S39_*`（报告）。

## 二、canonical 基因 manifest（D4 第一交付物）

`Mitoxyperilysis_Gene_Manifest_v1.0.csv`：80 基因 × 12 列（version / gene_symbol / hgnc_symbol / aliases / ensembl_gene_id / gene_name / module / module_size / arm / is_mt_gene / ARDS_vs_Control_log2FC / ARDS_vs_Control_padj）。

- **符号策略决策**（N1 遗留 2）：canonical 符号**维持现用常用名**（IP3R1、HSP60），不整体切换 HGNC——全部管线与结果表已用此写法，切换无分析收益且有断链风险；HGNC 官方名以 `hgnc_symbol` 列提供，别名以 `aliases` 列登记。
- **包内别名登记**：IP3R1↔ITPR1（gnomAD/AlphaMissense 层）、HSP60↔HSPD1（同前）、HSPA9↔GRP75（ATAC/转录本层）、CYCS↔CYTC（ATAC/TF 结合层）、SLC11A2↔DMT1（MR 层）。lint 符号校验按此归一化。
- **双臂归属**：上游塌陷臂 30（MAM 9+线粒体功能 10+铁死/铜死 11）、执行诱导臂 33（铁代谢 9+cell_death 13+氧化 11）、不在臂内 17（TF 11+自噬 6）——与 `singlecell_dissociation_test.py`（Table_S48，上游臂经别名归一化 30/30 全部命中，早期 28/30 系 IP3R1/HSP60 别名不一致）/ bulk Table_S5b / 主稿 §3 bulk log2FC 叙事为同一划分（**2026-08-15 N2 臂定义统一**：此前主稿叙事曾用"MAM+铁+线粒体+氧化 39 基因"旧口径、S5b 方法注曾按 Δz 符号事后归组把 TF/自噬并入上游组，均已统一为本 manifest `arm` 列；S5b 已增 `arm` 列）。
- 模块尺寸不变量：9/10/9/11/13/11/11/6 = 80（lint 强校验）。

## 三、版本戳（gene_set_version / score_version）

全部 173 个 canonical CSV 追加两常量列；21 个 TXT 报告追加版本页脚。取值：

| 版本号 | 含义 |
|---|---|
| `Mitoxy-80_v1.0` | 依赖 80 基因清单的分析（本 manifest 冻结版） |
| `PCD-6panel_v1.0` | S5 六种程序性死亡评分基因面板 |
| `score_genes_v1` | scanpy score_genes 逐细胞评分（S1、S6 Score_Summary、S7q、S48） |
| `ssGSEA_v1` | running-sum ssGSEA（β 参数见表内列；S5/S5b/S5c/S16f） |
| `bridge_v2` | Bridge 五维收敛（v2 = N1 别名归一化修复后重算） |
| `NA` | 不依赖基因集/评分版本（分组定义、全基因组 DEG 附属表、注释表等） |

**一致性校验**：lint 逐行核对每个 CSV 的版本列值 == manifest 记录值；TXT 页脚版本 == manifest。

## 四、登记的例外（lint 白名单，均附理由）

1. **TOMM70A**（层内扩展基因）：ATAC（S26f/S26j）、GWAS（S19c）、转录本（S41*）层的基因面板自带，非 80 清单成员——真实基因、真实数据，保留并登记。
2. **NADH**（伪符号）：GSE212865 微阵列探针注释伪影（R6 已登记）+ 转录本层沿袭（S41a–d/f）——保留源数据原样，不参与符号合法性判定。
3. **S40 对照 12 基因**：BAD/BAK1/BAX/BBC3/BCL2/BCL2L1/BCL2L11/BCL2L2/BID/MCL1/MCU/PMAIP1——正文 §1 Tier 比较的凋亡/MAM 钙对照，设计内扩展。
4. **S18 自选面板**：AlphaMissense 55 基因面板 = 80 清单成员 + 历史扩展基因（C2CD5、PRDX4、TOP2A 等），登记为 custom panel，不做 80 子集强校验。
5. **FDR/p 空值白名单**：S47b `FDR`（6 个恒定比例细胞类型相关未定义 → NaN，`FDR_note` 列逐行说明）；S8 `GSE212865_padj`/`meta_Fisher_padj`（基因不在该平台/某数据集 p 缺失，数据集缺席 NA）；S15a `fdr_q`（逐 IV 敏感度行，无 FDR 义务）。
6. **不加版本列的文件**：`GSE185263_groups.csv`（分组定义，读取方按列 merge）、`MR_bio_*.csv`（探索性 MR 结局侧，状态待定）、`GSE67530_beta_for_GrimAge.csv.gz`（输入数据）。

## 五、历史列名统一（R1 遗留①，数值零改动）

| 文件 | 改动 |
|---|---|
| `Table_S5b` | `ARDS_z`→`Sepsis_COVID_z`、`ARDS_pct_pos`→`Sepsis_COVID_pct_pos`、`delta`→`delta_Sepsis_COVID_vs_Control` |
| `Table_S5c` | `ARDS_z`→`Sepsis_COVID_z`、`delta_ARDS_Control`→`delta_Sepsis_COVID_vs_Control` |
| `Table_S9_Model_Performance.csv` | task 值 `ARDS_vs_Control`→`Sepsis_COVID_vs_Control`、`Sepsis_vs_ARDS`→`Sepsis_vs_Sepsis_COVID`；note 同步 |
| `Table_S9_Diagnostic_Report.txt` | Task A/B 标题与任务定位文字按 R1 正名（Sepsis_COVID 全血，非"ARDS 肺"） |
| `Bridge_Test_Result.csv` | `lung_log2FC`→`GSE185263_log2FC`（S46a 全血，历史"肺组织"误标）、`blood_log2FC`→`GSE212865_log2FC`（按数据集锚定） |

**保留不改的历史命名（登记备查）**：`Table_S46a` 的 `COVID_*` 列 = GSE185263 sepcv（Sepsis_COVID）——命名风格差异已在此登记，改动将级联大量读取脚本，收益为零。

**N3 溯源更正（2026-08-15，撤销本 manifest 旧登记）**：`Table_S2`/基因清单的 `ARDS_vs_Control_*` 列**不是** scRNA COVID_severe vs Healthy（本文件此前登记有误）。经数值指纹核对（ACTB≈13.0 为 bulk 计数尺度、EPCAM≈0.5 符合全血、样本设计 82/44 与 348/44 匹配 GSE185263、meanExpr=4.932 = log1p(452/3.3) 为 DESeq2 中位数比率归一化），该表实为 **GSE185263 全血 Bulk DESeq2 输出**（ARDS_vs_Control = 82 sepcv vs 44 ctrl；Sepsis_vs_Control = 348 sepsis vs 44 ctrl）。数值零改动，标签层面按 R1 先例更正；真·单细胞 DEG 由新增 `Table_S2b` 承担，旧全基因组 Fisher Meta 由新增 `Table_S8b` 替代（详见 `N3_Sample_Composition_Report.md`）。

## 六、lint 校验的不变量（计数与组定义）

GSE185263 = Sepsis 266 / Sepsis_COVID 82 / Control 44（共 392；groups、S16e、S5、S5_RAW 四处一致）；scRNA 合计 138,941 = GSE145926 83,952 + GSE158055 54,989（S1_Summary 与 S48 pooled 双处核验）；GSE145926 条件细胞数 26,138/7,050/50,764（N3 标签修复后，S2_QC 合计与 h5ad 一致）；S6 Scissor 组成比例逐组合计 100%；基因 manifest 80 行/8 模块/双臂 30-33-17；Bridge 80 行；**S2b 80 行（符号列严格校验）；S8b 80 行（符号列严格校验）**；S12b 23 TF；S13 12 KO；S15d 15 基因 34 cis+13 trans；S18 55 行；S46e 38 行；全部 p 值列 ∈[0,1]、FDR 列非空（白名单除外）；210 个 canonical 文件 SHA256 与本 manifest 一致。**（2026-08-16 S6 旧编号 PNG 归档后为 219 个。）**

## 七、未纳入本次冻结范围（显式登记）

1. **S46a/S46b/S46e 基于旧 38 基因清单生成**（N1 遗留 3）——是否按 80 清单重生成属 D1 数据决策，本次未改数值。
2. **`03_FIGURES/Figure*.csv`（140 个）**：历史出图用工作副本（独立 Figure 编号体系），**非补充表包成员**；其中 S16/S17/S18 等对应副本早于 2026-08-14 修复，出图前须自 canonical 表重新导出。
3. 叙述层重锚定（D1）、S40 gnomAD 提升（D2）等见质控报告路线图。

## 七a、旧版评分尺度标注（R5 处置，2026-08-15）

- `Table_S19h_CellType_Mitoxy_Score.csv` 与 `Table_S22b_scFOCAL_CellType_IC50.csv` 为**旧版评分尺度**（S19h：AUCell 式 0.46–0.90 值域，排序与现行 S1 score_genes 不一致；S22b：~20 量级评分 + 11 类含 Neu/Plasma/Mast 的细胞类型体系，与全稿 8 类 BALF 注释不兼容）——两表 `score_version` 列已标 `legacy_scale_v0_not_for_text`，**仅供历史参考，不用于正文数值引用**（R5"标注"选项；如需"重算"选项属后续工作量）。`Table_S1_rescore_report.txt` 标题已正名（实现为 score_genes，原"AUCell"名不符实）。

## 七b、登记的脚本侧同步（重跑一致性）

以下写入/读取脚本已同步至 canonical 命名（仅路径与标签文字，无算法改动）：`singlecell_dissociation_test.py`（S48）、`bridge_test_narrative.py`（GSE185263/GSE212865 列名+平台标签）、`build_Table_S5.py`（Sepsis_COVID 组名+新列名+叙述）、`NicheNet_Ligand_Receptor_Target_Analysis_v2.py`（S12b）、`Mitoxyperilysis_ARDS_GSE185263_WGCNA_v2.py` / `WGCNA_v2_rerun.py`（S49）、`annotate_cell_types_lung_atlas.py`（S26a）、`AlphaMissense_S18_Fix.py`（正式名输出）、`fix_s9_ml_model_v3_honest.py`（任务标签）、`rescore_mitoxy_80gene_aucell.py`（AUCell 注记）。

## 七c、Table_S36c 结构注释更正（R7 处置，2026-08-15）

`Table_S36c_Target_Protein_Info.csv` 的 Resolution/CoLigand/Function 列经 **RCSB 官方 REST API 逐条复核**后更正（备份 `archive/Table_S36c_pre_R7_backup.csv`，SHA 与冻结前原版一致）：

| 靶蛋白 | PDB | 旧（错误）值 | 新（RCSB 官方）值 |
|---|---|---|---|
| GPX4 | 5H5Q | 2.0Å；CoLigand "GXP" | **1.1Å**（X-ray）；CoLigand **GXpep-1**；Function 注明"结晶构建体含 Sec/Cys 替换" |
| VDAC1 | 2JK4 | 3.0Å | **4.1Å**（X-ray，Bayrhuber PNAS 2008）——注：质控报告 R7 曾称 2JK4 为"NMR 溶液结构"，经 RCSB 复核**该说法有误**（2JK4 为 X-ray 4.1Å；VDAC1 的 NMR 结构为 2K1N，未用于对接） |
| CANX | 1JHN | 未注物种/片段 | 2.9Å（不变）；Function 注明"**犬源**（*Canis lupus familiaris*）**腔域片段**，无跨膜/胞质段" |
| FTL | 2FHA | 2.5Å | **1.9Å**（X-ray） |
| TFRC | 1SUV | 3.2Å | **EM 7.5Å**（cryo-EM，人 TfR–transferrin 复合物；低分辨率模型，Function 已注明） |
| SOD2 | 1N0J | 2.2Å | 2.2Å（核验无误，未改） |
| SLC40A1 | AlphaFold | pLDDT=80.25 | 不变（预测结构标识原本正确） |

行数（7）与对接数值（S36/S37）零改动；仅注释列更正。主稿 Methods 分子对接段已同步披露结构证据分层（见文稿"深度学习与分子模拟"节）。manifest 重冻结 223 行 = 211 canonical + 12 archived，lint 全绿（50 项）。

## 八、重跑与再冻结

任何结果表变更后：`python build_D4_package.py`（仅对新文件补版本戳并重冻结 manifest）→ `python lint_package.py`（必须全绿）。冻结清单自身变更应升版本号（v1.1）并保留 v1.0 审计链。

## 九、附表编号迁移（2026-10-02，N7/N8 执行）

**动机**：最终投稿工作簿 `Supplementary_Tables_S1-S93.xlsx`（2026-10-01 版）把附表重排为与文稿 v5 一致的连续编号（S1–S93），本 manifest 的 `filename`/`s_number` 两列随之迁移（**原件 SHA256 不变，校验和无需改**）。

- 执行脚本：`N7_renumber_supp_tables_20261002.py`（重命名 + manifest/图数据溯源列同步）、`N8_register_wetlab_figdata_20261002.py`（湿实验图源与归档留档补登记）
- 逐文件旧→新对照：`02_SUPPLEMENTARY_TABLES/附表最终编号对照表_20261002.csv`；执行记录（含种子映射复核、备份与回滚）：`02_SUPPLEMENTARY_TABLES/RENUMBER_STATUS_20261002.md`
- 受影响：02 附表 **257 行** filename/s_number；**新增 canonical 21 行**（湿实验图源 `01_FIGURE_DATA_CSV/Main/Figure_8wA–8wH`、`Supplementary/Figure_S15A–G`、`Figure_S16A–F`）与 **archived 1 行**（`归档/Table_S26k_CellType_Scores.csv` 留档副本；canonical 件已改名 S79k）
- 组号示例：S5→S3、S6→S4、S7→S5、S15→S10、S18→S13、S26→S79、S48→S22、S49→S93、S55→S28（全量见上述 CSV）
- 7 个未纳入最终工作簿的文件保留原名（`Table_S3_Cell_Type_Composition_GSE158055/GSE171668_metadata`、`Table_S40_gnomAD_Constraint_Publication`、`Table_S5_PCD_Scores_RAW`、`Table_S7_cellphonedb_Raw`、`Table_S93/S93b_M17_*`）
- 本节之前各节的编号为**迁移前**口径（历史留痕，未逐一改写）
