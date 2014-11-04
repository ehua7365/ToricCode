"""
readData.py
"""
import sys
import math
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

def main():
    table = genfromtxt('data.csv', delimiter=',')
##    L = table[281,0]
##    N = table[281,3]
##    p = table[281:301,1].T
##    pl = table[281:301,4].T
    p = np.empty(0)
    pl = np.empty(0)
    for i in xrange(np.shape(table)[0]):
        if table[i,3] == 1000 and table[i,1] == .2:
            p = np.append(p,math.log(table[i,0]))
            pl = np.append(pl,math.log(table[i,4]))
    dp = (p[1]-p[0])/2
    #dpl = np.sqrt(pl*(1-pl)/N+(.5/N)**2)
    plt.errorbar(p,pl,yerr=0,xerr=dp,fmt=' ',c='0.4')#,label='$L = %i$'%L)
    c = np.polyfit(p,pl,6)
    print(c)
    P = np.poly1d(c)
    x = np.linspace(min(p),max(p),100)
    plt.plot(x,P(x),c='red')
    #plt.ylim(-0.05,1.05)
    plt.xlabel('Bit flip probability $p$')
    plt.ylabel('Unrecoverable probability $1-p_{recover}$')
    plt.title('Plot of $1-p_{recover}$ vs. $p$')
    plt.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    main()
