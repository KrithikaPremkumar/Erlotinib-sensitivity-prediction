from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr, spearmanr


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    return {"rmse": rmse, "mae": mae, "r2": r2}


def correlation_metrics(y_true, y_pred):
    if np.std(y_pred) < 1e-10:
        return {"pearson_r": float("nan"), "spearman_r": float("nan")}
    p_r, _ = pearsonr(y_true, y_pred)
    s_r, _ = spearmanr(y_true, y_pred)
    return {"pearson_r": p_r, "spearman_r": s_r}


def summarize_results(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    sort_cols = [c for c in ["rmse", "mae"] if c in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols).reset_index(drop=True)
    return df
