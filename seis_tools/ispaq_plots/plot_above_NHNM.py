#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

df = pd.read_csv("/home/conradb/git/ispaq/csv/pct_above_nhnm_NZ.CTZ.10.HHZ_2020-06-01_2021-01-31_PSDMetrics.csv")



df["start"] = pd.to_datetime(df["start"])
df = df.drop(df.columns[[0, 2, 3,]], axis=1)
print(df.head(20))

df["mean"] = df["value"].rolling(7, center=True).mean()
print(df.head(20))

df_melted = df.melt("start",var_name="value&mean",value_name="percent_value")

fig, ax = plt.subplots(figsize=(20,15))
#sns.lineplot(x="start", y="value" , data=df, marker="o")
sns.relplot(data=df_melted, x="start", y="percent_value", hue="value&mean",kind="line", height=4, aspect=.7)
ax.set(xlabel="Date", ylabel = "% above NHNM")
fig.autofmt_xdate()
plt.tight_layout()
plt.show()
