"""
readData.py
"""
import sys
import math
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

def main():
    table = genfromtxt('L20C.csv', delimiter=',')
    N = table[0,3]
    p = table[:,1].T
    pl = table[:,2].T
    dp = (p[1]-p[0])/2
    dpl = np.sqrt(pl*(1-pl)/N+(.5/N)**2)
    plt.errorbar(p,pl,yerr=dpl,xerr=dp,fmt=' ',c='0.4',label='$l = 20$')
    c = np.polyfit(p,pl,18)
    P = np.poly1d(c)
    x = np.linspace(min(p),max(p),100)
    plt.plot(x,P(x),c='red')
    plt.ylim(-0.05,1.05)
    plt.xlabel('Bit flip probability $p$')
    plt.ylabel('Unrecoverable probability $1-p_{recover}$')
    plt.title('Plot of $1-p_{recover}$ vs. $p$')
    plt.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    main()
