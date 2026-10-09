"""Work out the real hierarchy of the Rolle dataset.

Put this file next to reduced_dataset.pk, then run:
    pip install scipy
    python find_hierarchy.py
Copy everything it prints and paste it to Claude.

What it does:
1. Averages the 10-minute data to hourly (less noise, same structure).
2. Checks whether the aggregates add up to each other (e.g. all = S1 + S2).
3. For each aggregate, finds which of the 24 meters it is made of, using
   non-negative least squares: aggregate ~ w1*meter1 + ... + w24*meter24.
   If the aggregate is a plain sum, the weights come out as 1 (in) or 0 (out).
"""
import numpy as np
import pandas as pd
from scipy.optimize import nnls

df = pd.read_pickle("reduced_dataset.pk")
AGG = ["all", "S1", "S2", "S11", "S12", "S21", "S22"]
WEATHER = ["ghi_backwards", "temperature"]
meters = [c for c in df.columns if c not in AGG + WEATHER]
short = {m: f"m{i:02d}" for i, m in enumerate(meters, 1)}  # readable names

hourly = df[meters + AGG].astype(float).resample("1h").mean().dropna()
print(f"Hourly rows: {len(hourly)}  |  meters: {len(meters)}  |  aggregates: {AGG}")
print("\nMeter short names (used below):")
for m, s in short.items():
    print(f"  {s} = {m}")


def rel_gap(a, b):
    """Typical gap between two series, as % of the first."""
    return 100 * np.mean(np.abs(a - b)) / np.mean(np.abs(a))


print("\n=== 1. Do the aggregates add up to each other? (gap in %) ===")
checks = {
    "all vs S1+S2": (hourly["all"], hourly["S1"] + hourly["S2"]),
    "S1 vs S11+S12": (hourly["S1"], hourly["S11"] + hourly["S12"]),
    "S2 vs S21+S22": (hourly["S2"], hourly["S21"] + hourly["S22"]),
    "all vs S11+S12+S21+S22": (hourly["all"], hourly[["S11", "S12", "S21", "S22"]].sum(axis=1)),
    "all vs sum of 24 meters": (hourly["all"], hourly[meters].sum(axis=1)),
}
for name, (a, b) in checks.items():
    print(f"  {name:28s} gap = {rel_gap(a, b):6.2f}%")

print("\n=== 2. Which meters make up each aggregate? ===")
X = hourly[meters].to_numpy()
for agg in AGG:
    y = hourly[agg].to_numpy()
    w, _ = nnls(X, y)
    fitted = X @ w
    r2 = 1 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)
    used = [(short[m], round(wi, 2)) for m, wi in zip(meters, w) if wi > 0.05]
    print(f"\n  {agg}: R^2 = {r2:.4f}, gap = {rel_gap(y, fitted):.2f}%, meters used = {len(used)}")
    print("   ", ", ".join(f"{s}({wi})" for s, wi in used))

print("\nDone. Weights near 1.0 = that meter is part of the aggregate.")
