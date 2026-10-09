"""Make the clean hourly CSV for the pilot.

Put this file next to reduced_dataset.pk, then run:
    python make_hourly.py
It writes rolle_hourly.csv: one row per hour, columns in this order:
    time, all, S1, S2, S11, S12, S21, S22, m01 ... m24
(upper series first, then the 24 bottom meters, as bayesRecon expects).
"""
import pandas as pd

df = pd.read_pickle("reduced_dataset.pk")
UPPER = ["all", "S1", "S2", "S11", "S12", "S21", "S22"]
WEATHER = ["ghi_backwards", "temperature"]
meters = sorted(c for c in df.columns if c not in UPPER + WEATHER)
names = {m: f"m{i:02d}" for i, m in enumerate(meters, 1)}

hourly = (
    df[UPPER + meters].astype(float)
    .resample("1h").mean()
    .dropna()
    .rename(columns=names)
)
hourly.index.name = "time"
hourly.to_csv("rolle_hourly.csv")

# Save the meter-name key so we can always trace m01..m24 back to the dataset.
pd.Series(names, name="short_name").rename_axis("meter_id").to_csv("meter_names.csv")

print(f"Saved rolle_hourly.csv: {hourly.shape[0]} hours x {hourly.shape[1]} series")
print("From", hourly.index.min(), "to", hourly.index.max())
print(hourly.iloc[:3, :9].round(1))
