#!/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import linecache

plt.style.use("ggplot")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


data = np.loadtxt("data.txt", unpack = True)
params = linecache.getline("data.txt", 2)
params = params.split()[1:]

plt.plot(data[0], data[1], label = r'$\Delta_x(t)$')
plt.plot(data[0], data[2], label = r'$\Delta_y(t)$')
plt.plot(data[0], data[1][0]*np.ones(len(data[0])), label = r'$\Delta_0$')
plt.legend()
plt.xlim(np.min(data[0]), np.max(data[0]))
plt.xlabel(r'$t$')
plt.ylabel(r'$\Delta(t)$')
#plt.title(r'$\lambda={}$'.format(params[0]) + ",  " + r'\delta \lambda$={}$'.format(params[1]) + ",  "+ r'$\omega_D={}$'.format(params[2]) + ",  " + r'$\Delta_0={}$'.format(params[3]))
plt.savefig("statSol_{}.pdf".format(params))
plt.clf()








