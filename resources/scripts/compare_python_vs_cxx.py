import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
from matplotlib import colors 
import matplotlib.pyplot as plt

import numpy as np

data_cxx = np.loadtxt("../../examples/example_out1_cxx.txt", delimiter=",")
data_cxx_e = data_cxx[:,0]
data_cxx_r = data_cxx[:,1]

data_py = np.loadtxt("../../src/python/example_out1_python.txt", delimiter=",")
data_py_e = data_py[:,0]
data_py_r = data_py[:,1]

plt.figure()
fig, ax = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)

ax[0].plot(data_cxx_e, data_cxx_r, "b-", label="cxx version")
ax[0].plot(data_py_e, data_py_r, "r-", label="python version")
ax[0].set_ylabel("transmittance")
ax[0].set_title("E^-2 transmittance at 180 deg")
ax[0].grid()
ax[0].legend(loc="best")

ax[1].plot(data_cxx_e, data_cxx_r/data_cxx_r, "b-", label="cxx/cxx")
ax[1].plot(data_py_e, data_py_r/data_cxx_r, "r-", label="python/cxx ")
ax[1].set_xlabel("log10(primary energy)")
ax[1].set_ylabel("ratio")
ax[1].set_ylim((0.8, 1.4))
ax[1].grid()
ax[1].legend(loc="best")

plt.savefig("python_vs_cxx.png")


