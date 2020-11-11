# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sinh, tanh, cosh, abs, log
from matplotlib.ticker import ScalarFormatter
from scipy import integrate
from scipy.special import ellipk, ellipe
import os
import csv

tau_array = np.arange(0.1,5.0,0.1)
tau_c = 2.0/(log(1+sqrt(2)))
print("Critical temperature : ", tau_c)

############### parameter ################
n = 100
Nspin = n*n
# 0:熱浴法 1:メトロポリス法
METHOD = 0
############### parameter ################

###############################
if(METHOD == 0):
  DATA_PATH = "./data/heat_bath/spin%d"%n
  SAVE_PATH = "./figs/heat_bath/spin%d"%n
if(METHOD == 1):
  DATA_PATH = "./data/metropolis/spin%d"%n
  SAVE_PATH = "./figs/metropolis/spin%d"%n
df_result = pd.read_csv("%s/result.csv"%DATA_PATH)
###############################

def my_K(k):
  ret = []
  for k_ in k:
    f = lambda theta : 1.0/np.sqrt(1.0-k_*k_*np.sin(theta)**2)
    ret.append(integrate.quad(f, 0.0, pi/2.0)[0])
  return np.array(ret)

def my_E(k):
  ret = []
  for k_ in k:
    f = lambda theta : np.sqrt(1.0-k_*k_*np.sin(theta)**2)
    ret.append(integrate.quad(f, 0.0, pi/2.0)[0])
  return np.array(ret)

def set_ax(ax):
  # ax.set_title("Voight Profile", fontsize=24)
  # ax.set_ylabel("Voight Profile", fontsize=21)
  # ax.set_ylabel(r"$\phi$(x)", fontsize=21)
  # ax.set_xlabel("x", fontsize=21)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  # ax.tick_params(labelsize=21)
  # ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
  # ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
  # ax.yaxis.offsetText.set_fontsize(20)

def coth(x):
  return 1/tanh(x)

def u_ana(tau):
  k = 2.0/(cosh(2.0/tau)*coth(2.0/tau))
  k2 = 2.0*tanh(2.0/tau)*tanh(2.0/tau)-1
  # tmp = 1+2.0/pi*k2*ellipk(k)
  tmp = 1+2.0/pi*k2*my_K(k)
  return -1.0/tau*coth(2.0/tau)*tmp

def c_ana(tau):
  k = 2.0/(cosh(2.0/tau)*coth(2.0/tau))
  k2 = 2.0*tanh(2.0/tau)*tanh(2.0/tau)-1
  # tmp = 2.0*ellipk(k)-2.0*ellipe(k)-(1-k2)*(pi/2.0+k2*ellipk(k))
  tmp = 2.0*my_K(k)-2.0*my_E(k)-(1-k2)*(pi/2.0+k2*my_K(k))
  return 2.0/(tau*tau*pi)*coth(2.0/tau)*coth(2.0/tau)*tmp

def m_ana(tau):
  ret = [] 
  for t in tau:
    if(t < tau_c):
      tmp = 1-sinh(2.0/t)**(-4.0)
      ret.append(tmp**(1.0/8.0))
    else:
      ret.append(0)
  return ret

def u_plot():
  tau = df_result["T"]
  u = df_result["e"]
  fig = plt.figure(figsize=(8, 8))
  fig.subplots_adjust(left=0.2)
  ax = fig.add_subplot(111)
  set_ax(ax)
  ax.axvline(x=tau_c,color="r",ls="--",label=r"$\tau_c$")
  ax.set_xlabel(r"$\tau$",fontsize=20)
  ax.set_ylabel(r"u",fontsize=20)

  ax.plot(tau_array, -2/tau_array, "--", label=r"2/$\tau$")
  ax.plot(tau_array, u_ana(tau_array), label="analytic")
  ax.plot(tau, u, ".-", label="numeric")
  ax.legend()
  plt.savefig("%s/u.png"%SAVE_PATH)
  # plt.show()


def c_plot():
  tau = df_result["T"]
  c = df_result["c"]

  fig = plt.figure(figsize=(8, 8))
  fig.subplots_adjust(left=0.2)
  ax = fig.add_subplot(111)
  set_ax(ax)
  ax.axvline(x=tau_c,color="r",ls="--",label=r"$\tau_c$")
  ax.set_xlabel(r"$\tau$",fontsize=20)
  ax.set_ylabel(r"c",fontsize=20)
  ax.plot(tau_array, c_ana(tau_array), label="analytic")
  ax.plot(tau, c, ".-", label="numeric")
  ax.legend()
  plt.savefig("%s/c.png"%SAVE_PATH)
  # plt.show()

def m_plot():
  tau = df_result["T"]
  m = df_result["m"]

  fig = plt.figure(figsize=(8, 8))
  fig.subplots_adjust(left=0.2)
  ax = fig.add_subplot(111)
  set_ax(ax)
  ax.axvline(x=tau_c,color="r",ls="--",label=r"$\tau_c$")
  ax.set_xlabel(r"$\tau$",fontsize=20)
  ax.set_ylabel(r"m",fontsize=20)
  ax.plot(tau_array, m_ana(tau_array), label="analytic")
  ax.plot(tau, m, ".-" ,label="numeric")
  ax.legend()
  plt.savefig("%s/m.png"%SAVE_PATH)
  # plt.show()

def K_E_compare():
  df = pd.read_csv("./data/ellips.csv")
  k = df["k"]
  K_pp = df["K"]
  E_pp = df["E"]
  fig = plt.figure(figsize=(8, 8))
  # fig.subplots_adjust(left=0.2)
  ax = fig.add_subplot(111)
  set_ax(ax)
  ax.set_xlabel(r"k",fontsize=20)
  ax.set_ylabel("K(k),E(k)",fontsize=20)
  ax.plot(k, ellipk(k), "b", label="scipy K(k)")
  ax.plot(k, my_K(k), "b--" ,label="my K(k)")
  ax.plot(k, K_pp, "b.-" ,label="C++ K(k)")
  ax.plot(k, ellipe(k), "r", label="scipy E(k)")
  ax.plot(k, my_E(k), "r--" ,label="my E(k)")
  ax.plot(k, E_pp, "r.-" ,label="C++ E(k)")
  ax.legend()
  # plt.show()
  plt.savefig("./figs/EK_compare.png")

u_plot()
c_plot()
m_plot()
# K_E_compare()
plt.show()