# M13 RCD 去混淆审计报告

- 日期：2026-08-26；依据：M10_M13_M14_pre_registration_20260824.md §2（osf.io/ETVMJ）
- 签名目录：MSigDB v7.5.1 HALLMARK/C2.CP(KEGG+REACTOME)/C5.GO 关键词筛选 + 10<=n<=600 + 去重 -> 249 个签名
- 臂定义（去循环化主口径）：外部本体 上游=GOCC 线粒体被膜/内膜/基质/RESPIRASOME+GOBP OXPHOS（1105 基因）；执行=GOBP 焦亡/炎症小体装配+GOCC 炎症小体复合物+REACTOME 焦亡+gasdermin/CASP1 家族（71 基因）
- 髓系基因：Monaco 29 型 argmax 于髓系类型（11999 基因）

## 1. 预注册回归（主检验）


- 模型：签名疾病效应 g（GSE185263 Sepsis_COVID vs Control）~ 执行臂占比 + 髓系基因占比
- **beta_exec=+0.719（se=0.490），单侧 p=0.07169**
- beta_myeloid=+2.309（p=6.27e-11）；R^2=0.206，n=249
- LOSO：beta_exec 方向一致 249/249（100%，阈值>=90%）
- 零模型（B=10,000 表达量匹配随机等大小基因集）：null beta_exec=+0.609；bootstrap 经验 p=0.4196
- 零模型（myeloid 臂，2026-08-29 溯源补录）：null beta_myeloid=**1.612±0.063**（bootstrap 2000 次重拟合，seed=0，复现自 `M13_step2_null_bootstrap.py`，输出落盘 `03_LOGS/M13_null_bootstrap_myeloid_log.txt`）；myeloid bootstrap 经验 p=**0.0005**（观测 beta_myeloid=2.309；B=10,000 次重拟合稳健核验：0/10000 达 2.309 → p=0.0001，零分布 1.613±0.062，数值稳定）
- 敏感性（manifest 臂定义）：beta_exec=+2.321

## 2. 反向审计


| 论文 | 策展集 | overlap/可测 | 期望 | Fisher p | 匹配零模型 p |
|---|---|---|---|---|---|
| PMC13189447_melanoma_BAX_BAK1_BID | BAX,BAK1,BID | 2/3 | 1.78 | 0.639 | 0.9098 |
| PMC13189447_melanoma_+mTORC2 | BAX,BAK1,BID,MTOR,RICTOR | 3/5 | 2.97 | 0.672 | 0.9413 |
| PMC13189447_melanoma_+mTORC2_full | BAX,BAK1,BID,MTOR,RICTOR,MLST8,MAPKAP1 | 4/7 | 4.16 | 0.698 | 0.9516 |

（措辞纪律：只陈述基因集一致性，不评价他人论文结论对错。）

## 3. 判定门 M13


- 判定门：预注册回归显著 + LOSO 稳定（>=90% 方向一致）+ 外部本体定义下同向。判定见正文。