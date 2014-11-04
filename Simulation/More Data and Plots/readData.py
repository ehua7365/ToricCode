"""
readData.py
"""
import sys
import math
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import datetime

def main():
    table = genfromtxt('data2.csv', delimiter=',')
    Lmax = int(max(table[:,0]))
    colormap = plt.cm.gist_ncar
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0,0.9,Lmax*2)])
    n = np.empty(0)
    totalTime = np.empty(0)
    for L in xrange(2,Lmax):
        x = np.empty(0)
        y = np.empty(0)
        N = np.empty(0,dtype=int)
        t = np.empty(0)
        for row in xrange(table.shape[0]):
            if int(table[row,0]) == L:
                x = np.append(x,table[row,1])
                y = np.append(y,table[row,2])
                N = np.append(N,int(table[row,3]))
                t = np.append(t,table[row,4])
        pmax = max(x)
        pstep = pmax/len(x)
        z = np.polyfit(x,y,4)
        p = np.poly1d(z)
        xp = np.linspace(0,pmax,100)
        plt.plot(xp,p(xp))
        plt.errorbar(x,y,xerr=pstep/2,yerr=np.sqrt(y*(1-y)/N),fmt='.',label='$L = %i$'%L)
        n = np.append(n,L*L)
        totalTime = np.append(totalTime,sum(t))
    plt.ylim(-0.05,1.05)
    plt.xlabel('Bit flip probability $p$')
    plt.ylabel('Unrecoverable probability $1-p_{recover}$')
    plt.title('Plot of $1-p_{recover}$ vs. $p$')
    plt.legend(loc='upper left')
    fileName = 'Figure%s.pdf'%datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    plt.savefig(fileName)
    plt.show()

    plt.plot(n,totalTime)
    plt.xlabel('Number of Qubits $n = L^2$')
    plt.ylabel('Computation Time $t_C$ (s)')
    plt.title('Plot of $t_C$ vs. $n$')
    fileName = 'CompTime%s.pdf'%datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    plt.savefig(fileName)
    plt.show()

if __name__ == '__main__':
    main()
