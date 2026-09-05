# P13 待办 1：OSF 三个注册链接人工验证记录（2026-09-04）

> 对应 `111提质修复_20260902.md` 待办 1。
> 验证方式（双重）：① 作者（用户）于 2026-09-04 在浏览器逐一打开三个注册页，确认公开状态、
> 页面标题与 DOI，并将结果回传记录在案；② 同日经 OSF 官方 API 独立佐证
> （`curl --ssl-no-revoke https://api.osf.io/v2/registrations/<ID>/`，本机 WebFetch 被拦但 curl 可达），
> title / date_registered / public / identifiers 四字段与人工核验一致。

## 结论：三项全部通过

| 注册 | OSF 页面标题（用户回传，逐字） | DOI | 对应模块 | 核对结果 |
|---|---|---|---|---|
| C7RYD | Frozen Analysis Plan v1.0 — Phase 0 Recomputation (ARDS Mitoxyperilysis-aligned program) | 10.17605/OSF.IO/C7RYD | P0 冻结分析计划 | 公开 ✓ 标题对应 ✓ DOI 逐字 ✓ |
| ETVMJ | Secondary-data preregistration of the composition-hub adjudication (M10), regulated-cell-death deconfounding audit (M13), and perturbation triangulation (M14) for the two-arm dissociation in critical-illness transcriptomes (Mitoxyperilysis ARDS project) | 10.17605/OSF.IO/ETVMJ | M10/M13/M14 | 公开 ✓ 标题对应 ✓ DOI 逐字 ✓ |
| 98CM3 | Preregistration of M11/M12: temporal sequencing and intervention response of the composition-independent mitochondrial-infrastructure suppression axis | 10.17605/OSF.IO/98CM3 | M11/M12 时间序与干预响应 | 公开 ✓ 标题对应 ✓ DOI 逐字 ✓ |

## 与待办 1 核对清单的逐项对应

1. **公开访问（不跳登录墙）**：三个链接用户均确认公开可访问。
2. **标题与模块对应**：三标题与预期模块逐一相符（见上表）。
3. **date_registered（98CM3 应为 2026-08-27）**：**已由 OSF API 独立确认**（2026-09-04）：

   | 注册 | API date_registered (UTC) | 与既有留痕的一致性 |
   |---|---|---|
   | C7RYD | 2026-08-20T05:00:21Z | 与 P0 重算冻结时间线一致（osf.io/3e92j 旧版之后的正式注册版） |
   | ETVMJ | 2026-08-24T14:12:51Z | 与 ETVMJ 注册文档日期（2026-08-24）一致 |
   | 98CM3 | 2026-08-27T05:45:15Z | **逐秒吻合**注册文档自述 2026-08-27T05:45:15Z；UTC+8 即 13:45，与表单审计留痕的确认邮件时间戳一致；早于全部 M11/M12 结局统计量计算 |
4. **DOI 与两稿 Data availability 一致**：`111文稿_v5_英文版.md`（L607）与
   `归档/111文稿_v5.md`（L597）均为 10.17605/OSF.IO/C7RYD、/ETVMJ、/98CM3（大写），
   与 OSF 页面逐字一致；两稿 Table 4 溯源节（英文 L133 / 中文 L129）同。

## 大小写订正说明（98cm3 → 98CM3）

OSF 规范写法为 **98CM3**（小写 URL 可达，但项目内统一按大写）。

- **已订正（2026-09-04）**：
  - `111提质修复_20260902.md` 待办 1 链接列表（原文误写小写）；
  - `02_SUPPLEMENTARY_TABLES/README_附表索引.md` L133（活文档）。
- **不订正（冻结/历史留痕，保留原样）**：
  - `M11_M12_pre_registration_20260827_draft.md` L4——已注册冻结件，内容与 OSF
    存档一致优先，不回改本地副本；
  - `04_AUDIT_GOVERNANCE/M11M12_Analysis_Report_20260827.md`、
    `05_RELEASE_GITHUB/repo/governance/` 同名镜像——带日期审计报告按惯例原样保留；
  - `03_LOGS/M11M12_OSF_form_audit_20260827.md`、`M11M12_OSF_form_field_corrections_20260827.md`——
    当日操作日志；
  - `_intermediate/DOI_replace_log.txt` 及 `DOI_*.py`——历史替换操作留痕；
  - `M11M12_step4_export.py`（含 `05_RELEASE_GITHUB/repo/analysis_scripts/` 镜像）L160——
    冻结导出脚本；若待办 4 发布前重跑该脚本，建议顺手把该行字符串改为大写。

## 其他

- 截图未归档（用户口头确认）。如投稿前需要，可补充浏览器截图存入本目录。
- 结论：待办 1 关闭；两稿 Data availability 与 Table 4 脚注无需任何修改。
