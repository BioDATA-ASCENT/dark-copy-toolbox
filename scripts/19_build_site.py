#!/usr/bin/env python3
"""S9 — assemble the GitHub Pages site in docs/:
  docs/index.html    the Dark Copy Toolbox (reports/explorer/index.html, built by 14_export_explorer.py)
  docs/outline.html  the project outline (explorer/outline.html)
  docs/img/*.png     figures the outline uses (from reports/)
"""
import os
import re
import shutil

from cfg import REPO_ROOT
from common import status

DOCS = os.path.join(REPO_ROOT, "docs")
os.makedirs(os.path.join(DOCS, "img"), exist_ok=True)
shutil.copy(os.path.join(REPO_ROOT, "reports", "explorer", "index.html"), os.path.join(DOCS, "index.html"))
outline = open(os.path.join(REPO_ROOT, "explorer", "outline.html")).read()
open(os.path.join(DOCS, "outline.html"), "w").write(outline)
for img in sorted(set(re.findall(r'src="img/([^"]+)"', outline))):
    shutil.copy(os.path.join(REPO_ROOT, "reports", img), os.path.join(DOCS, "img", img))
open(os.path.join(DOCS, ".nojekyll"), "w").close()
status("S9_site", "DONE", "docs/index.html + docs/outline.html")
