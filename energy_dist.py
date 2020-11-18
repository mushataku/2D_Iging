import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os




FILE_PATH = 'energy_histo'

fig = plt.figure(figsize=(14, 8))
# fig.subplots_adjust(left=0.80)
ax = fig.add_subplot(111)

T = ["0.100","0.400", "0.800", "2.000", "5.000"]
labels = [r"$\tau$ = 0.1", r"$\tau$ = 0.4",
          r"$\tau$ = 0.8", r"$\tau$ = 2.0", r"$\tau$ = 5.0"]
Energys = []

for t in T:
  print(t)
  df = pd.read_csv("data/spin100/energy_dist/spin100_T" + t + ".csv")
  Energy = df["Energy"]
  Energys.append(Energy)

ax.hist(Energys, bins=33, ec="black", align="left", range=(-100,32), 
stacked=False, label=labels, rwidth=30)


ax.grid(linestyle="dotted")
ax.xaxis.set_tick_params(direction='in')
ax.yaxis.set_tick_params(direction='in')
ax.tick_params(labelsize=21)  # 軸目盛の数字のサイズ
ax.legend(fontsize=21)
ax.set_yscale("log")
plt.savefig("fig/"+FILE_PATH+".png")
plt.savefig("report_tex/"+FILE_PATH+".png")
plt.show()
