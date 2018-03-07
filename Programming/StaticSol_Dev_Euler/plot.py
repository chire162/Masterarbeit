import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.txt", unpack = True)

plt.plot(data[0], data[1])
plt.plot(data[0], data[2])
plt.savefig("plot.pdf")
plt.show()