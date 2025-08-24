
"""
atp_analysis.py

Reusable ATP analysis toolkit for luminometer runs with time-drifting Tris blanks.

Key features
------------
1) Fits the ATP standard curve FIRST (on raw integrals from standards),
   then models Tris drift and applies it as a time-varying intercept.
2) Clear comments, input validation, and predictable behaviors.
3) Packaged as a class for import and use in Jupyter or scripts.
4) Optional: read a CSV with per-sample weights/volumes to back-calculate original concentrations.
5) Optional: merge with a metadata file (CSV/TSV).
6) Optional: specify blank names and compute k·SD threshold above blanks.
7) Optional: plot against a value/column from the metadata.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple, Union, Callable
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import time, datetime


def _to_time(obj: Union[str, time, pd.Timestamp]) -> time:
    """Convert to naive datetime.time."""
    if isinstance(obj, time):
        return obj
    if isinstance(obj, pd.Timestamp):
        return obj.to_pydatetime().time()
    if isinstance(obj, str):
        # Try common formats; fallback to pandas
        for fmt in ("%m/%d/%y %H:%M", "%Y-%m-%d %H:%M:%S", "%H:%M:%S", "%H:%M:%S.%f"):
            try:
                return datetime.strptime(obj, fmt).time()
            except ValueError:
                continue
        return pd.to_datetime(obj).to_pydatetime().time()
    raise TypeError(f"Cannot convert type {type(obj)} to time")

def _ensure_parent_dir(path: str) -> None:
    d = os.path.dirname(path) or "."
    os.makedirs(d, exist_ok=True)

def _seconds_since(start: time, t: time) -> float:
    return (t.hour - start.hour) * 3600 + (t.minute - start.minute) * 60 + (t.second - start.second)


def _extract_standard_concentration(sample_name: str) -> Optional[float]:
    """Parse sample/index name for a leading numeric standard concentration."""
    if not sample_name:
        return None
    s = str(sample_name).strip()
    m = re.match(r"^(?:\d+|\d*\.\d+)", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


@dataclass
class FitResult:
    slope: float
    intercept: float

    def predict(self, x: np.ndarray, through_origin: bool = False) -> np.ndarray:
        """Predict y for plotting the standard curve."""
        x = np.asarray(x, dtype=float)
        if through_origin:
            return self.slope * x
        return self.slope * x + self.intercept


class ATPAnalyzer:
    """
    End-to-end ATP analysis.

    Typical flow:
        analyzer = (ATPAnalyzer(data_csv, time_csv)
                      .integrate()
                      .fit_standard_curve()
                      .fit_tris_drift(separate=True, split_seconds=1000)
                      .apply_corrections_and_quantify(extract_vol=4.0, sample_vol=weights_df, sample_unit="g")
                      .compute_blank_threshold(blank_names=['Blank','Blank.1'], k=3.0)
                      .merge_metadata(meta_path, right_key='#SampleID'))
        analyzer.plot_standard_curve('ATP_Standard_Curve.png', through_origin=True)
        analyzer.plot_vs_metadata(x='Age', y='avg_concentration', save_path='atp_vs_age.png')
        analyzer.save_outputs(prefix='outputs/atp_')
    """

    def __init__(self, data_csv: str, time_csv: str) -> None:
        self.data_csv = data_csv
        self.time_csv = time_csv

        self.integrals_df: Optional[pd.DataFrame] = None
        self.standard_fit: Optional[FitResult] = None
        self.tris_drift_func: Optional[Callable[[np.ndarray], np.ndarray]] = None
        self._tris_t0: Optional[time] = None

        self.samples_wide: Optional[pd.DataFrame] = None
        self.grouped_samples: Optional[pd.DataFrame] = None
        self.merged_meta: Optional[pd.DataFrame] = None

        self._standards_grouped_for_plot: Optional[pd.DataFrame] = None

    # ------------------------
    # Data integration loader
    # ------------------------
    def integrate(self) -> "ATPAnalyzer":
        """
        Load time series and compute trapezoidal integrals for each sample column.
        The timestamps CSV must contain a single row; its headers must match sample names.
        """
        df = pd.read_csv(self.data_csv)
        if "Time" not in df.columns:
            raise ValueError("data_csv must include a 'Time' column representing elapsed time (seconds).")
        # Drop accidental non-sample columns like 'Unnamed: XX' and columns that are entirely NaN (except 'Time')
        drop_cols = [c for c in df.columns if c != "Time" and (str(c).startswith("Unnamed") or df[c].isna().all())]
        if drop_cols:
            df = df.drop(columns=drop_cols)

        time_df = pd.read_csv(self.time_csv)
        if time_df.shape[0] != 1:
            time_df = time_df.iloc[[0]]
        time_map: Dict[str, time] = {}
        for k, v in time_df.iloc[0].to_dict().items():
            try:
                time_map[k] = _to_time(v)
            except Exception:
                time_map[k] = np.nan

        integrals: Dict[str, float] = {}
        for col in df.columns:
            if col == "Time":
                continue
            integrals[col] = np.trapz(df[col].values, df["Time"].values)

        out = pd.DataFrame.from_dict(integrals, orient="index", columns=["Integral"])
        out["Time"] = out.index.map(time_map)
        self.integrals_df = out
        return self

    # ------------------------
    # Fitting steps
    # ------------------------
    def _split_standards_and_samples(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if self.integrals_df is None:
            raise RuntimeError("Call .integrate() first.")
        df = self.integrals_df.copy()
        df["Standard_Concentration"] = df.index.to_series().apply(_extract_standard_concentration)
        standards = df[df["Standard_Concentration"].notna()].copy()
        nonstandards = df[df["Standard_Concentration"].isna()].copy()
        return standards, nonstandards

    def fit_standard_curve(self) -> "ATPAnalyzer":
        """Fit linear standard curve on RAW integrals (pre-Tris-correction)."""
        standards, _ = self._split_standards_and_samples()
        if standards.empty:
            raise ValueError("No standards detected. Ensure sample/index names begin with a number (e.g., '100', '0.5').")
        grouped = (
            standards.groupby("Standard_Concentration", as_index=False)["Integral"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={"mean": "mean_integral", "std": "sd_integral", "count": "n"})
        )
        x = grouped["Standard_Concentration"].to_numpy(float)
        y = grouped["mean_integral"].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        self.standard_fit = FitResult(float(slope), float(intercept))
        self._standards_grouped_for_plot = grouped.copy()
        return self

    def fit_tris_drift(self, separate: bool = False, split_seconds: Optional[float] = None) -> "ATPAnalyzer":
        """Fit a Tris drift model (single line or two segments)."""
        if self.integrals_df is None:
            raise RuntimeError("Call .integrate() first.")
        df = self.integrals_df.copy()
        tris_df = df[df.index.str.contains("Tris", case=False, na=False)].copy()
        if tris_df.empty:
            # zero-drift fallback
            self.tris_drift_func = lambda t: np.zeros_like(np.asarray(t, dtype=float))
            # set t0 for later
            valid = df["Time"].dropna().tolist()
            if not valid:
                raise ValueError("No valid timestamps found to apply drift correction.")
            self._tris_t0 = min(valid, key=lambda x: (x.hour, x.minute, x.second))
            return self

        valid_times = tris_df["Time"].dropna().tolist()
        if not valid_times:
            raise ValueError("Tris rows found but no valid timestamps were parsed.")
        t0 = min(valid_times, key=lambda x: (x.hour, x.minute, x.second))
        tris_sec = np.array([_seconds_since(t0, _to_time(tt)) for tt in tris_df["Time"]], dtype=float)
        tris_int = tris_df["Integral"].to_numpy(float)

        if not separate:
            a, b = np.polyfit(tris_sec, tris_int, 1)
            self.tris_drift_func = lambda t: a * np.asarray(t, dtype=float) + b
            self._tris_model = {"mode": "single", "a": float(a), "b": float(b)}
        else:
            if split_seconds is None:
                raise ValueError("split_seconds must be provided when separate=True.")
            m1 = tris_sec <= split_seconds
            m2 = ~m1
            if m1.sum() < 2 or m2.sum() < 2:
                raise ValueError("Not enough Tris points on one side of split_seconds to fit two lines.")
            a1, b1 = np.polyfit(tris_sec[m1], tris_int[m1], 1)
            a2, b2 = np.polyfit(tris_sec[m2], tris_int[m2], 1)

            def piecewise(t: Union[float, np.ndarray]) -> np.ndarray:
                tt = np.asarray(t, dtype=float)
                out = np.empty_like(tt)
                mask1 = tt <= split_seconds
                mask2 = tt > split_seconds
                out[mask1] = a1 * tt[mask1] + b1
                out[mask2] = a2 * tt[mask2] + b2
                return out

            self.tris_drift_func = piecewise
            self._tris_model = {
                "mode": "piecewise",
                "a1": float(a1), "b1": float(b1),
                "a2": float(a2), "b2": float(b2),
                "split": float(split_seconds),
            }

        self._tris_t0 = t0
        return self

    # --------------------------------------------
    # Quantification, blanks, metadata, and plots
    # --------------------------------------------
    def apply_corrections_and_quantify(
        self,
        extract_vol: float,
        sample_vol: Optional[Union[float, pd.DataFrame]] = None,
        sample_volume_key: str = "Sample_Volume",
        base_sample_key: str = "Base_Sample",
        sample_unit: Optional[str] = None,   # e.g., "g" or "mL" for labeling
    ) -> "ATPAnalyzer":
        """
        Subtract predicted Tris baseline at each sample's timestamp,
        then convert corrected luminescence to concentration using the standard slope.
        """
        if self.integrals_df is None:
            raise RuntimeError("Call .integrate() first.")
        if self.standard_fit is None:
            raise RuntimeError("Call .fit_standard_curve() first.")
        if self.tris_drift_func is None:
            self.tris_drift_func = lambda t: np.zeros_like(np.asarray(t, dtype=float))
            valid = self.integrals_df["Time"].dropna().tolist()
            if not valid:
                raise ValueError("No valid timestamps available to apply drift correction.")
            self._tris_t0 = min(valid, key=lambda x: (x.hour, x.minute, x.second))

        # seconds since start
        tsec = self.integrals_df["Time"].apply(lambda tt: _seconds_since(self._tris_t0, _to_time(tt)) if pd.notna(tt) else np.nan)
        drift_pred = self.tris_drift_func(tsec.to_numpy(dtype=float))

        df = self.integrals_df.copy()
        df["Predicted_Tris_Luminescence"] = drift_pred
        df["Corrected_Luminescence"] = (df["Integral"] - df["Predicted_Tris_Luminescence"]).clip(lower=0.0)
        self.integrals_df = df
        # quantify
        _, nonstandards = self._split_standards_and_samples()
        samples = nonstandards[~nonstandards.index.str.contains("Tris", case=False, na=False)].copy()
        samples["Base_Sample"] = samples.index.to_series().astype(str).str.split(".", n=1).str[0]

        # Use ONLY slope (Tris acts as dynamic intercept)
        m = self.standard_fit.slope
        samples["Concentration_ng_per_mL"] = samples["Corrected_Luminescence"] / m

        # Back-calc to original sample
        if sample_vol is None:
            samples["Total_ng_in_extract"] = samples["Concentration_ng_per_mL"] * float(extract_vol)
        elif isinstance(sample_vol, (int, float)):
            unit = sample_unit or "sample_unit"
            colname = f"Original_Sample_Concentration_ng_per_{unit}"
            samples[colname] = samples["Concentration_ng_per_mL"] * float(extract_vol) / float(sample_vol)
        elif isinstance(sample_vol, pd.DataFrame):
            if not {"Base_Sample", sample_volume_key}.issubset(sample_vol.columns):
                raise ValueError(f"sample_vol DataFrame must contain ['Base_Sample', '{sample_volume_key}']")
            unit = sample_unit or "sample_unit"
            vol_map = sample_vol.set_index("Base_Sample")[sample_volume_key]
            samples["Sample_Volume"] = samples["Base_Sample"].map(vol_map)
            colname = f"Original_Sample_Concentration_ng_per_{unit}"
            samples[colname] = samples["Concentration_ng_per_mL"] * float(extract_vol) / samples["Sample_Volume"].astype(float)
        else:
            raise TypeError("sample_vol must be None, a float, or a pandas DataFrame.")

        self.samples_wide = samples

        # Group to mean/stdev per Base_Sample
        conc_candidates = [c for c in samples.columns if c.startswith("Original_Sample_Concentration_ng_per_")]
        if conc_candidates:
            conc_col = conc_candidates[0]
        else:
            # fall back to extract concentration
            conc_col = "Concentration_ng_per_mL"
        grouped = (
            samples.groupby("Base_Sample", as_index=False).agg(
                avg_concentration=(conc_col, "mean"),
                std_concentration=(conc_col, "std"),
                n=("Corrected_Luminescence", "count"),
            )
        )
        self.grouped_samples = grouped
        return self

    def compute_blank_threshold(self, blank_names: Iterable[str], k: float = 3.0) -> "ATPAnalyzer":
        """Compute blank mean/SD and add above_k_sd_blanks to grouped table."""
        if self.samples_wide is None or self.grouped_samples is None:
            raise RuntimeError("Call .apply_corrections_and_quantify() first.")
        blanks = self.samples_wide[self.samples_wide.index.isin(blank_names)].copy()
        if blanks.empty:
            blanks = self.samples_wide[self.samples_wide["Base_Sample"].isin(blank_names)].copy()
        if blanks.empty:
            self.blank_mean = np.nan
            self.blank_sd = np.nan
            self.blank_threshold = np.nan
            self.grouped_samples["above_k_sd_blanks"] = np.nan
            return self

        conc_candidates = [c for c in blanks.columns if c.startswith("Original_Sample_Concentration_ng_per_")]
        if conc_candidates:
            conc_col = conc_candidates[0]
        elif "Concentration_ng_per_mL" in blanks.columns:
            conc_col = "Concentration_ng_per_mL"
        else:
            raise ValueError("No concentration column found in blanks.")

        self.blank_mean = float(blanks[conc_col].mean())
        self.blank_sd = float(blanks[conc_col].std(ddof=1))
        self.blank_threshold = self.blank_mean + k * self.blank_sd
        self.grouped_samples["above_k_sd_blanks"] = self.grouped_samples["avg_concentration"] > self.blank_threshold
        return self

    def merge_metadata(
        self,
        meta: Union[str, pd.DataFrame],
        left_key: str = "Base_Sample",
        right_key: str = "#SampleID",
        strip_trailing_period_in_meta_ids: bool = True,
    ) -> "ATPAnalyzer":
        """Merge grouped per-sample stats with metadata file or DataFrame."""
        if self.grouped_samples is None:
            raise RuntimeError("Call .apply_corrections_and_quantify() first.")
        if isinstance(meta, str):
            if meta.endswith((".tsv", ".tab", ".txt")):
                meta_df = pd.read_csv(meta, sep="\t")
            else:
                meta_df = pd.read_csv(meta)
        else:
            meta_df = meta.copy()
        if strip_trailing_period_in_meta_ids and right_key in meta_df.columns:
            meta_df[right_key] = meta_df[right_key].astype(str).str.rstrip(".")
        merged = pd.merge(self.grouped_samples.copy(), meta_df, left_on=left_key, right_on=right_key, how="left")
        self.merged_meta = merged
        return self

    # -------------
    # Plot helpers
    # -------------
    def plot_standard_curve(self, save_path: Optional[str] = None, through_origin: bool = False) -> "ATPAnalyzer":
        """Plot standards (mean±SD) with fitted line; optionally show through-origin line."""
        if self.standard_fit is None or self._standards_grouped_for_plot is None:
            raise RuntimeError("Call .fit_standard_curve() first.")
        grouped = self._standards_grouped_for_plot
        x = grouped["Standard_Concentration"].to_numpy(float)
        y = grouped["mean_integral"].to_numpy(float)
        yerr = grouped["sd_integral"].to_numpy(float)
        xx = np.linspace(float(x.min()), float(x.max()), 200)
        yy = self.standard_fit.predict(xx, through_origin=through_origin)

        plt.figure(figsize=(8, 6))
        plt.errorbar(x, y, yerr=yerr, fmt='o', capsize=5, label="Standards (mean ± SD)")
        label = "Linear fit (through origin)" if through_origin else f"Linear fit: y = {self.standard_fit.slope:.3g}x + {self.standard_fit.intercept:.3g}"
        plt.plot(xx, yy, label=label)
        plt.xlabel("[ATP] (ng/mL)")
        plt.ylabel("Raw integral")
        plt.title("ATP Standard Curve (fit before Tris correction)")
        plt.legend()
        if save_path:
            _ensure_parent_dir(save_path)
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return self

    def plot_vs_metadata(
        self,
        x: str,
        y: str = "avg_concentration",
        yerr: str = "std_concentration",
        save_path: Optional[str] = None,
        add_blank_threshold: bool = True,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
    ) -> "ATPAnalyzer":
        """Scatter with error bars using merged metadata."""
        if self.merged_meta is None:
            raise RuntimeError("Call .merge_metadata() first.")
        plot_df = self.merged_meta.dropna(subset=[x, y]).copy()
        plt.figure(figsize=(10, 7))
        plt.errorbar(plot_df[x], plot_df[y], yerr=plot_df.get(yerr, None), fmt="o", ecolor="lightgray", elinewidth=2, capsize=4)
        if add_blank_threshold and hasattr(self, "blank_threshold") and pd.notna(getattr(self, "blank_threshold", np.nan)):
            plt.axhline(self.blank_threshold, linestyle="--", label="Blank threshold", linewidth=1)
        plt.xlabel(xlabel or x)
        plt.ylabel(ylabel or y)
        plt.title(f"{y} vs {x}")
        plt.legend()
        if save_path:
            _ensure_parent_dir(save_path)
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return self
    def plot_tris_drift(self, save_path: Optional[str] = None) -> "ATPAnalyzer":
        """
        Plot Tris integrals vs seconds since start, with the fitted drift lines.
        - If single-line model: draw one line.
        - If piecewise model: draw TWO separate line segments split at `split`.
        """
        if self.integrals_df is None or self._tris_t0 is None or getattr(self, "_tris_model", None) is None:
            raise RuntimeError("Call .fit_tris_drift() first.")

        tris_df = self.integrals_df[self.integrals_df.index.str.contains("Tris", case=False, na=False)].copy()
        if tris_df.empty:
            raise ValueError("No Tris samples found.")

        # X,Y data
        tris_df["tsec"] = tris_df["Time"].apply(lambda tt: _seconds_since(self._tris_t0, _to_time(tt)) if pd.notna(tt) else np.nan)
        x = tris_df["tsec"].to_numpy(dtype=float)
        y = tris_df["Integral"].to_numpy(dtype=float)

        # Build plot
        plt.figure(figsize=(8, 6))
        plt.scatter(x, y, marker="o", label="Observed Tris", zorder=3)

        model = self._tris_model
        if model["mode"] == "single":
            a, b = model["a"], model["b"]
            xx = np.linspace(np.nanmin(x), np.nanmax(x), 200)
            yy = a * xx + b
            plt.plot(xx, yy, label=f"Drift fit: y={a:.3g}·t + {b:.3g}", zorder=2)
        else:
            a1, b1, a2, b2, split = model["a1"], model["b1"], model["a2"], model["b2"], model["split"]
            # left segment
            xx1 = np.linspace(np.nanmin(x), split, 100)
            yy1 = a1 * xx1 + b1
            # right segment
            xx2 = np.linspace(split, np.nanmax(x), 100)
            yy2 = a2 * xx2 + b2
            plt.plot(xx1, yy1, label=f"Drift fit (≤ {split:.0f}s): y={a1:.3g}·t + {b1:.3g}", zorder=2)
            plt.plot(xx2, yy2, label=f"Drift fit (> {split:.0f}s): y={a2:.3g}·t + {b2:.3g}", zorder=2)
            # vertical split guide
            plt.axvline(split, linestyle="--", linewidth=1, label="split", zorder=1)

        plt.xlabel("Time since start (s)")
        plt.ylabel("Raw integral")
        plt.title("Tris baseline drift")
        plt.legend()

        if save_path:
            _ensure_parent_dir(save_path)
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return self
    # -------------
    # I/O helpers
    # -------------
 
   

    def save_outputs(
        self,
        prefix: str = "atp_",
        save_integrals: bool = True,
        save_samples_wide: bool = True,
        save_grouped: bool = True,
        save_merged_meta: bool = True,
    ) -> "ATPAnalyzer":
        """Save CSV outputs if present."""
        # ensure directory exists
        out_dir = os.path.dirname(prefix) or "."
        os.makedirs(out_dir, exist_ok=True)

        if save_integrals and self.integrals_df is not None:
            self.integrals_df.to_csv(f"{prefix}integrals.csv")
        if save_samples_wide and self.samples_wide is not None:
            self.samples_wide.to_csv(f"{prefix}samples_wide.csv", index=True)
        if save_grouped and self.grouped_samples is not None:
            self.grouped_samples.to_csv(f"{prefix}grouped.csv", index=False)
        if save_merged_meta and self.merged_meta is not None:
            self.merged_meta.to_csv(f"{prefix}merged_meta.csv", index=False)
        return self