# -*- coding: utf-8 -*-
"""Solve NCBI PMC cloudpmc-viewer PoW and download DataSheet1.docx for PMC13376262."""
import hashlib
import http.cookiejar
import re
import sys
import urllib.request

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"

cj = http.cookiejar.CookieJar()
opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))


def get(url):
    return opener.open(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=120)


file_url = "https://pmc.ncbi.nlm.nih.gov/articles/instance/13376262/bin/DataSheet1.docx"
page = get(file_url).read().decode(errors="replace")
print("page len:", len(page))
m = re.search(r'POW_CHALLENGE\s*=\s*"([^"]+)"', page)
if not m:
    print(page[-800:])
    sys.exit("no challenge found")
challenge = m.group(1)
difficulty = int(re.search(r'POW_DIFFICULTY\s*=\s*"([^"]+)"', page).group(1))
name = re.search(r'POW_COOKIE_NAME\s*=\s*"([^"]+)"', page).group(1)
print("challenge len:", len(challenge), "difficulty:", difficulty, "cookie:", name)

prefix = "0" * difficulty
nonce = 0
while True:
    h = hashlib.sha256((challenge + str(nonce)).encode()).hexdigest()
    if h.startswith(prefix):
        break
    nonce += 1
print("nonce:", nonce, "hash:", h)

cj.set_cookie(http.cookiejar.Cookie(
    version=0, name=name, value=f"{challenge},{nonce}", port=None, port_specified=False,
    domain="pmc.ncbi.nlm.nih.gov", domain_specified=True, domain_initial_dot=False,
    path="/", path_specified=True, secure=True, expires=None, discard=False,
    comment=None, comment_url=None, rest={}, rfc2109=False))

data = get(file_url).read()
print("bytes:", len(data), "magic:", data[:8])
assert data[:2] == b"PK", "not a docx zip"
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\PMC13376262_DataSheet1.docx"
open(out, "wb").write(data)
print("saved ->", out)

