# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, pi
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin, cosh, abs
from matplotlib.ticker import ScalarFormatter
from scipy import integrate
from scipy.special import *
import os
import csv

SHOW_2 = 0

spinnum = "100"

df_E = pd.read_csv("data/spin"+spinnum+"/energy.csv")

T = df_E["T"]
E_numeric = df_E["numerical"]
E_analytic = df_E["analytic"]
E_analytic2 = df_E["analytic2"]
error = df_E["error"]
error2 = df_E["error2"]

E_FILE_PATH = "spin"+spinnum+"energy"
e_FILE_PATH = "spin"+spinnum+"e_error"

#plt.rcParams["font.family"] = "Times New Roman"

fig_E = plt.figure(figsize=(12, 8))
# fig.subplots_adjust(left=0.80)
Eax = fig_E.add_subplot(111)

Eax.set_title("N = "+spinnum, fontsize=26)
Eax.plot(T, E_numeric,"r-o", label="numeric")
Eax.plot(T, E_analytic,"b", label="analytic")
if(SHOW_2 == 1):
  Eax.plot(T, E_analytic2,"--", label="analytic2")
Eax.set_ylabel(r"energy e", fontsize=24)
Eax.set_xlabel(r"Temperature $\tau$", fontsize=24)
Eax.grid(linestyle = "dotted")
Eax.xaxis.set_tick_params(direction='in')
Eax.yaxis.set_tick_params(direction='in')
Eax.set_xticks(np.arange(0, 5.1, 0.5))
Eax.tick_params(labelsize=20)  # 軸目盛の数字のサイズ
Eax.legend(fontsize=20)
# os.mkdir("fig")
plt.savefig("fig/"+E_FILE_PATH+".png")
if(SHOW_2 == 0):
  plt.savefig("report_tex/"+E_FILE_PATH+".png")
else:
  plt.savefig("report_tex/"+E_FILE_PATH+"2.png")

fig_e = plt.figure(figsize=(12, 8))
# fig.subplots_adjust(left=0.80)
eax = fig_e.add_subplot(111)

eax.set_title("N = "+spinnum, fontsize=26)
eax.plot(T, error,"r-o", label="error")
if(SHOW_2 == 1):
  eax.plot(T, error2,"b-o", label="error2")
eax.set_ylabel(r"Rerative Error [%]", fontsize=24)
eax.set_xlabel(r"Temperature $\tau$", fontsize=24)
eax.grid(linestyle = "dotted")
eax.xaxis.set_tick_params(direction='in')
eax.yaxis.set_tick_params(direction='in')
eax.set_xticks(np.arange(0, 5.1, 0.5))
eax.set_yscale("log")
eax.tick_params(labelsize=20)  # 軸目盛の数字のサイズ
eax.legend(fontsize=20)
# os.mkdir("fig")
plt.savefig("fig/"+e_FILE_PATH+".png")
if(SHOW_2 == 0):
  plt.savefig("report_tex/"+e_FILE_PATH+".png")
else:
  plt.savefig("report_tex/"+e_FILE_PATH+"2.png")

#plt.show()
