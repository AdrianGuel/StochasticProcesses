import sdepy
from sdepy import *  # safe and handy for interactive sessions
import numpy as np
import scipy
import matplotlib.pyplot as plt  # optional, if plots are needed
from numpy import sqrt

@integrate
def my_process(t, x, u=1., k=10., sigma=1):
    return {'dt': -k*(x-u), 'dw': sigma}

#myp = kfunc(my_process)
timeline = np.linspace(0., 1, 101)
np.random.seed(5)  # make doctests predictable

x = my_process(x0=-1, vshape=1,paths=100000, steps=1)(timeline)
for k in timeline:
    y = x(k)[0]  # 0-th component of x at time t=1
    a = montecarlo(y, bins=50)
    ygrid = np.linspace(y.min(), y.max(), 200)
    plt.plot(ygrid, a.pdf(ygrid))
#gr = plt.plot(ygrid, a.pdf(ygrid, method='interp', kind='nearest'))
plt.ylim(top=5)
plt.show()  # doctest: +SKIP


