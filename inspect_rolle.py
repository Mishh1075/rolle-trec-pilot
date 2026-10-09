"""Inspect the Rolle distribution-grid dataset (Zenodo 4549296).

Put this file in the same folder as reduced_dataset.pk, then run:
    python inspect_rolle.py
Copy everything it prints and paste it to Claude.
"""
import os
import sys

import pandas as pd

FILE = sys.argv[1] if len(sys.argv) > 1 else "reduced_dataset.pk"

if not os.path.exists(FILE):
    sys.exit(f"Cannot find {FILE} in {os.getcwd()}. Put this script next to the file.")

obj = pd.read_pickle(FILE)
print("TYPE:", type(obj))


def describe(df, name="data"):
    print(f"\n=== {name} ===")
    print("shape (rows, columns):", df.shape)
    print("index:", type(df.index).__name__, df.index.min(), "->", df.index.max())
    if isinstance(df.index, pd.DatetimeIndex) and len(df.index) > 1:
        print("time step:", pd.Series(df.index).diff().mode()[0])
    if isinstance(df.columns, pd.MultiIndex):
        print("column levels:", df.columns.names)
        for i in range(df.columns.nlevels):
            vals = list(df.columns.get_level_values(i).unique())
            print(f" level {i} ({len(vals)} values):", vals)
    else:
        print(f"columns ({len(df.columns)}):", list(df.columns))
    print("missing values:", int(df.isna().sum().sum()))
    with pd.option_context("display.max_columns", 8, "display.width", 120):
        print(df.head(3))


def walk(o, name="data", depth=0):
    if isinstance(o, pd.DataFrame):
        describe(o, name)
    elif isinstance(o, pd.Series):
        describe(o.to_frame(), name)
    elif isinstance(o, dict) and depth < 3:
        print(f"\n{name}: dict with keys {list(o.keys())}")
        for k, v in o.items():
            walk(v, f"{name}[{k!r}]", depth + 1)
    else:
        print(f"\n{name}:", type(o), str(o)[:300])


walk(obj)
