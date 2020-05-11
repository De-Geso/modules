import numpy as np
from numpy import loadtxt, zeros, linspace
import matplotlib
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import itertools


path = '../../data/interpolation/'

x, y = loadtxt(path+'interpolation.dat', unpack=True)
xpol, ypol, dypol = loadtxt(path+'polint.dat', unpack=True)
xrat, yrat, dyrat = loadtxt(path+'ratint.dat', unpack=True)


fig1, ax1 = plt.subplots()
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Polynomial Interpolation')

fig2, ax2 = plt.subplots()
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('Rational Interpolation')


color = next(ax1._get_lines.prop_cycler)['color']
ax1.scatter(x, y, color=color)
ax1.plot(xpol, ypol, color=color)
ax1.fill_between(xpol, ypol-dypol, ypol+dypol, color=color, alpha=0.5)

ax2.scatter(x, y, color=color)
ax2.plot(xrat, yrat, color=color)
ax2.fill_between(xrat, yrat-dyrat, yrat+dyrat, color=color, alpha=0.5)


fig1.tight_layout()
fig2.tight_layout()


fig1.savefig('polint.pdf')
fig2.savefig('ratint.pdf')


plt.show()
