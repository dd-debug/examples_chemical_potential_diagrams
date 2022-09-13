'''
Created on Nov 15, 2021

@author: jiadongc
'''
import matplotlib.pyplot as plt
import numpy as np
import math
# create 1000 equally spaced points between -10 and 10
x = np.linspace(0, 7, 1000)
# print(x)
# calculate the y value for each element of the x vector
y = 4*x**2 + x-50


fig, ax = plt.subplots()
fig1,ax1 = plt.subplots()

# fig = plt.figure(figsize=(9.2, 7))
# ax = plt.gca()
# fig1 = plt.figure(figsize=(9.2, 7))
# ax1 = plt.gca()

ax.plot(x, y, linewidth=3, zorder=10)
ax.set_xlim([0,5])
ax.set_ylim([-100,10])
point = []
primalpoints = []
for i,c in [[3,"lightgray"],[2.5,"lightgray"],[1,"r"],[2,"g"]]:
    xx = i
    k = 8*xx + 1
    yy = 4*xx**2 + xx-50
    b=yy-k*xx
    z = k*x+b
    print(f'z={k}x {b}')
    point.append([k,b,c])
    if c == "lightgray":
        ax.plot(x, z, linewidth=1.5,c=c,zorder=5)
    else:
        ax.plot(x, z, linewidth=2,c=c,zorder=8)
        ax.scatter(xx,yy,c=c,s=80,zorder=12)
        primalpoints.append([xx,yy,c])
# plt.show()

pointsx = []
pointsy = []

for i in x:
    if i<5:
        xx = i
        k = 8*xx + 1
        yy = 4*xx**2 + xx-50
        b=yy-k*xx
        z = k*x+b
        pointsx.append(k)
        pointsy.append(b)
ax1.plot(pointsx, pointsy,linewidth=3,zorder=2)
for i in point:
    ax1.scatter(i[0],i[1],c=i[2],s=80,zorder=3)

a = np.linspace(0, 30, 10000)
for xx,yy,c in primalpoints:
    b=-(xx*a-yy)
    ax1.plot(a,b,c=c,linewidth=2,zorder=1)  
ax1.set_xlim([0,30])
ax1.set_ylim([-100,10])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
ax1.tick_params(axis='y', length = 2)
ax.tick_params(axis='y', length = 2)
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
plt.show()