#!/usr/bin/python3
# -*- coding: utf-8 -*-

import pandas as pd
import csv
import numpy as np
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
import seaborn as sns
from plotnine import *
from plotnine.data import * 
import scipy
from scipy import interpolate
import numpy
from numpy import linspace
import math
from math import sin
import re

df = pd.read_csv('itog.csv', sep = ',') # в df забиваем таблицу


x1 = df.iloc[:,[1]] #выбираем столбцы и обозначаем x,y
y1 = df.iloc [:,[2]]
x2 = np.array(x1)
y2 = np.array(y1)
print (x1, y1)

#a, b = np.polyfit(x2, y2, deg=1)
#y_est = a * x2 + b
#y_err = x2.std() * np.sqrt(1/len(x2) +
                         # (x2 - x2.mean())**2 / np.sum((x2 - x2.mean())**2))

#fig, ax = plt.subplots()
#ax.plot(x2, y_est, '-')
#ax.fill_between(x2, y_est - y_err, y_est + y_err, alpha=0.2)
#ax.plot(x2, y2, 'o', color='tab:brown')


#df.plot(kind='scatter', x='Steps', y='Mean', color = 'g')
#df.vlines(b1)
#plt.title("Frequency distribution during asexual reproduction", fontsize=18)
#plt.grid(True)
#plt.show()