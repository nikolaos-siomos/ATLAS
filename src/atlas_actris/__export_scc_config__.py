#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 18:20:24 2025

@author: nikos
"""

import csv
from pathlib import Path
from typing import Dict, Any, Optional

def export_lidar_config(data: Dict[str, Dict[Any, Dict[str, Any]]],
                        path: str,
                        checksum: Optional[str] = None) -> None:
    """
    Export the nested dictionary (from parse_lidar_config) back to a file
    with the original sectioned CSV structure.
    
    :param data: dict in form { section_name: {row_key: {col: val, ...}, ...}, ... }
    :param path: output filename
    :param checksum: optional checksum string (without "CHECKSUM" keyword)
    """
    with open(Path(path), "w", encoding="utf-8", newline="") as f:
        for idx, (section_name, rows) in enumerate(data.items()):
            # write section name
            f.write(section_name + "\n")

            if rows and isinstance(rows, dict):
                # figure out header fields from the first row
                sample_row = next(iter(rows.values()))
                fieldnames = list(sample_row.keys())
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in rows.values():
                    writer.writerow({k: _stringify(v) for k, v in row.items()})
            # blank line after each section
            f.write("\n")

        if checksum:
            f.write(f"CHECKSUM {checksum}\n")

def _stringify(v):
    """Convert Python types back to the textual CSV form."""
    if v is True: 
        return "True"
    if v is False: 
        return "False"
    if v is None: 
        return ""
    return str(v)

# -----------------------------
# Example CLI usage
# -----------------------------
if __name__ == "__main__":
    import sys, json
    from lidar_parser import parse_lidar_config  # reuse parser from before

    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.csv> <output.csv> [checksum]")
        sys.exit(1)

    infile, outfile = sys.argv[1], sys.argv[2]
    checksum = sys.argv[3] if len(sys.argv) > 3 else None

    cfg = parse_lidar_config(infile)
    export_lidar_config(cfg, outfile, checksum)
    print(f"Exported {len(cfg)} sections to {outfile}")