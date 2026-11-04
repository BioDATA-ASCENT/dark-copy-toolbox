#!/usr/bin/env python3
"""Read config/config.yaml.

CLI:    python scripts/cfg.py urls.grch38_fasta      -> prints the value
Module: from cfg import CFG, path                    -> dict + absolute-path helper
"""
import os
import sys

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
with open(os.path.join(REPO_ROOT, "config", "config.yaml")) as fh:
    CFG = yaml.safe_load(fh)


def path(*parts):
    """Absolute path under data_dir (created if it is a directory-like path)."""
    return os.path.join(REPO_ROOT, CFG["data_dir"], *parts)


def logs(*parts):
    return os.path.join(REPO_ROOT, CFG["logs_dir"], *parts)


def get(dotted):
    value = CFG
    for key in dotted.split("."):
        value = value[key]
    return value


if __name__ == "__main__":
    v = get(sys.argv[1])
    if isinstance(v, (list, dict)):
        print(yaml.safe_dump(v, default_flow_style=True).strip())
    else:
        print(v)
