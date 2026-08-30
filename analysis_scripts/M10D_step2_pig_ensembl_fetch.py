# -*- coding: utf-8 -*-
"""
M10D_step2_pig_ensembl_fetch.py — 猪同源 Ensembl 敏感性复算·缓存回填（人→sus_scrofa）
======================================================================================
背景：M10D 主口径执行时（2026-08-27）Ensembl REST 批量端点 404/单符号端点持续 503，
     仅有 24 个符号缓存（23 有效 + 1 个 503 残留），one2one 可用 8 条 → 猪层降级为
     GPL16524 注释符号逐字口径（51/63）。本脚本在服务恢复窗口重试回填缓存，
     供 M10D_step3_pig_ensembl_sensitivity.py 以 Ensembl 一对一映射复算 D3 敏感性。

口径（与主脚本一致，先于数值冻结）：
  - 端点：GET /homology/symbol/homo_sapiens/{SYM}?target_species=sus_scrofa;type=orthologues
    （批量 POST 端点本次实测 500，服务端错误 → 不用批量）。
  - 缓存格式与主脚本读取逻辑兼容：{"data":[{"homologies":[...]}]}；"data" 为空的
    合法响应（如 MT-* 线粒体基因无 Compara 同源组）同样落盘。
  - 重试：每符号 ≤5 次、读超时 300 s、指数退避；4 线程；断点续传（已有有效缓存跳过）。
  - 日志：03_LOGS/M10D_pig_ensembl_fetch_log.txt（先冻结口径再请求）。
"""
import json
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import mdi_lib as L

ROOT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(ROOT, "_intermediate", "m10d_pig_homology")
LOGP = os.path.join(ROOT, "03_LOGS", "M10D_pig_ensembl_fetch_log.txt")
MAX_TRY = 5
TIMEOUT = 300
WORKERS = 4

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def flush_log():
    with open(LOGP, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")

def cache_valid(path):
    try:
        r = json.loads(open(path, encoding="utf-8", errors="replace").read())
        return isinstance(r, dict) and isinstance(r.get("data"), list)
    except Exception:
        return False

def fetch_one(sym):
    url = ("https://rest.ensembl.org/homology/symbol/homo_sapiens/%s"
           "?target_species=sus_scrofa;type=orthologues;content-type=application/json" % sym)
    last_err = "?"
    for attempt in range(1, MAX_TRY + 1):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "Mozilla/5.0 (scientific analysis)"})
            raw = urllib.request.urlopen(req, timeout=TIMEOUT).read()
            r = json.loads(raw)
            if not (isinstance(r, dict) and isinstance(r.get("data"), list)):
                raise ValueError("unexpected payload")
            with open(os.path.join(CACHE, sym + ".json"), "w", encoding="utf-8") as f:
                json.dump(r, f)
            n = len(r["data"])
            return (sym, "ok", n, None)
        except Exception as e:
            last_err = str(e)[:120]
            if attempt < MAX_TRY:
                time.sleep(5 * attempt)
    return (sym, "FAIL", None, last_err)

def main():
    note("== M10D 猪同源 Ensembl 缓存回填（敏感性复算前置）2026-08-27")
    note("口径冻结：GET 单符号端点；批量 POST 本次实测 500（服务端）不采用；"
         "重试≤%d、超时%d s、%d 线程、断点续传；缓存格式与 M10D 主脚本读取逻辑一致；"
         "data=[] 的合法空响应（如 MT-*）同样落盘。" % (MAX_TRY, TIMEOUT, WORKERS))
    man, arms, _ = L.load_manifest()
    arm_all = sorted(set(arms["upstream_collapse"] + arms["execution_induction"]))
    todo = [s for s in arm_all if not cache_valid(os.path.join(CACHE, s + ".json"))]
    note("臂基因 %d 个；已有有效缓存 %d；本次待取 %d 个：%s"
         % (len(arm_all), len(arm_all) - len(todo), len(todo), ", ".join(todo)))
    t0 = time.time()
    ok, empty, fail = [], [], []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(fetch_one, s): s for s in todo}
        for fut in as_completed(futs):
            sym, status, n, err = fut.result()
            if status == "ok":
                (ok if n > 0 else empty).append(sym)
                note("[ok %d/%d] %s (homologies=%s) %.0fs"
                     % (len(ok) + len(empty), len(todo), sym, n, time.time() - t0))
            else:
                fail.append(sym)
                note("[FAIL] %s (%s)" % (sym, err))
    note("完成：ok=%d（%s）| 合法空=%d（%s）| 失败=%d（%s）| 耗时 %.0f s"
         % (len(ok), ", ".join(ok), len(empty), ", ".join(empty),
            len(fail), ", ".join(fail), time.time() - t0))
    flush_log()

if __name__ == "__main__":
    main()
