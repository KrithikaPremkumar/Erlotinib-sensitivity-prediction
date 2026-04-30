from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_ccle_gct(gct_path: Path) -> pd.DataFrame:
    """Load a GCT file into a tidy dataframe indexed by gene symbol."""
    raw = pd.read_csv(
        gct_path,
        sep="\t",
        skiprows=2,
        low_memory=False,
        dtype={"Name": str, "Description": str},
    )
    if "Description" not in raw.columns:
        raise ValueError("Expected GCT columns to include 'Description'.")

    # Keep the first occurrence of gene symbol for stable mapping.
    raw = raw.dropna(subset=["Description"]).drop_duplicates(subset=["Description"], keep="first")
    raw = raw.set_index("Description")
    sample_cols = [c for c in raw.columns if c != "Name"]
    expr = raw[sample_cols].apply(pd.to_numeric, errors="coerce")
    return expr


def clean_cellline_name(name: str) -> str:
    if not isinstance(name, str):
        return ""
    cleaned = name.strip().upper().replace("-", "").replace("_", "").replace(" ", "")
    return cleaned


def infer_ccle_sample_mapping(expr_df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for col in expr_df.columns:
        if "_" in col:
            cell_line = col.split("_")[0]
            tissue = col.rsplit("_", maxsplit=1)[-1]
        else:
            cell_line = col
            tissue = "UNKNOWN"
        records.append(
            {
                "ccle_sample_col": col,
                "cell_line_raw": cell_line,
                "cell_line_clean": clean_cellline_name(cell_line),
                "tissue_hint": tissue,
            }
        )
    return pd.DataFrame(records)


def load_gdsc_response(
    gdsc_path: Path,
    drug_name: str = "Erlotinib",
    ic50_column_hint: Optional[str] = None,
) -> pd.DataFrame:
    suffix = gdsc_path.suffix.lower()
    if suffix == ".csv":
        df = pd.read_csv(gdsc_path)
    elif suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(gdsc_path)
    else:
        raise ValueError(f"Unsupported GDSC file format: {gdsc_path.suffix}")
    lower_cols = {c.lower(): c for c in df.columns}

    drug_col = None
    for c in df.columns:
        if "drug" in c.lower() and ("name" in c.lower() or "id" not in c.lower()):
            drug_col = c
            break
    if drug_col is None:
        raise ValueError("Could not infer drug-name column in GDSC file.")

    ic50_col = ic50_column_hint
    if ic50_col is None:
        # GDSC2 often has both IC50 and LN_IC50; prefer LN_IC50 (already log µM).
        normalized = {c: c.lower().replace(" ", "") for c in df.columns}
        for c, lc in normalized.items():
            if "ln_ic50" in lc:
                ic50_col = c
                break
        if ic50_col is None:
            for c, lc in normalized.items():
                if "ic50" in lc:
                    ic50_col = c
                    break
    if ic50_col is None:
        raise ValueError("Could not infer IC50 column in GDSC file.")

    cell_line_col = None
    for candidate in df.columns:
        lc = candidate.lower()
        if "cell" in lc and "line" in lc:
            cell_line_col = candidate
            break
    if cell_line_col is None and "cosmic_id" in lower_cols:
        cell_line_col = lower_cols["cosmic_id"]
    if cell_line_col is None:
        raise ValueError("Could not infer cell line identifier column.")

    filtered = df[df[drug_col].astype(str).str.upper() == drug_name.upper()].copy()
    filtered = filtered[[cell_line_col, ic50_col]].dropna()
    filtered.columns = ["cell_line_raw", "ic50_raw"]
    filtered["cell_line_clean"] = filtered["cell_line_raw"].astype(str).map(clean_cellline_name)
    filtered["ic50_raw"] = pd.to_numeric(filtered["ic50_raw"], errors="coerce")
    filtered = filtered.dropna(subset=["ic50_raw", "cell_line_clean"])
    return filtered


def align_ccle_with_gdsc(
    expr_df: pd.DataFrame,
    gdsc_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    mapping = infer_ccle_sample_mapping(expr_df)
    # Avoid column collision: both sides have cell_line_raw; merge would produce _x/_y.
    gdsc_for_merge = gdsc_df.rename(columns={"cell_line_raw": "gdsc_cell_line_raw"})
    merged = mapping.merge(gdsc_for_merge, on="cell_line_clean", how="inner")
    if merged.empty:
        raise ValueError("No matched cell lines found between CCLE and GDSC after cleaning.")

    selected_cols = merged["ccle_sample_col"].tolist()
    x = expr_df[selected_cols].T
    # Index = CCLE cell-line token (matches prefix of ccle_sample_col).
    x.index = merged["cell_line_raw"].values
    y = pd.Series(merged["ic50_raw"].values, index=x.index, name="ic50")
    return x, y, merged


def safe_log_ic50(y: pd.Series) -> pd.Series:
    """Map IC50 to a modeling target. Log only for strictly positive concentration IC50.

    GDSC often ships **LN_IC50** (natural log µM), which can be negative; do not log again.
    """
    y = pd.to_numeric(y, errors="coerce")
    valid = y.dropna()
    if valid.empty:
        return y
    if (valid > 0).all():
        min_positive = float(valid.min())
        offset = min_positive / 10 if min_positive > 0 else 1e-6
        return np.log(y + offset)
    return y.copy()
