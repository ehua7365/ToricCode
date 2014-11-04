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
    table = genfromtxt('data.csv', delimiter=',')
    Lmax = int(max(table[:,0]))
    colormap = plt.cm.gist_ncar
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0,0.9,Lmax*2)])
    n = np.empty(0)
    totalTime = np.empty(0)
    for L in xrange(2,Lmax):
        x = []
        y = []
        N = []
        t = []
        for row in xrange(table.shape[0]):
            if int(table[row,0]) == L:
                x.append(table[row,1])
                y.append(table[row,2])
                N.append(int(table[row,3]))
                t.append(table[row,4])
        # Merge duplicate x-values
        index = 0
        while index < len(x)-1:
            for i in xrange(index+1,len(x)-2):
                if index < i and i < len(x):
                    if x[index] == x[i]:
                        x.pop(i)
                        y[index] = (y[index]*N[index]+y.pop(i)*N[i])/(N[index]+N[i])
                        t[index] = (t[index]/N[i]+t.pop(i)/N[index])*(N[index]+N[i])
                        N[index] += N.pop(i)
            index += 1
        # Convert to NumPy arrays
        x = np.array(x)
        y = np.array(y)
        N = np.array(N)
        t = np.array(t)
        # Calculate polynomial fit
        pmax = max(x)
        pstep = pmax/len(x)
        z = np.polyfit(x,y,4)
        p = np.poly1d(z)
        xp = np.linspace(0,pmax,100)
        # Plot fit and points
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

    z = np.polyfit(n,totalTime,10)
    print(z)
    p = np.poly1d(z)
    x = np.linspace(min(n),max(n),100)

    plt.plot(x,p(x),c='red')
    plt.errorbar(n,totalTime,fmt='+',c='blue')
    plt.xlabel('Number of Qubits $n = L^2$')
    plt.ylabel('Computation Time $t_C$ (s)')
    plt.title('Plot of $t_C$ vs. $n$')
    fileName = 'CompTime%s.pdf'%datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    plt.savefig(fileName)
    plt.show()

    lnn = np.log(n)
    lntt = np.log(totalTime)
    z = np.polyfit(lnn[15:],lntt[15:],1)
    print(z)
    p = np.poly1d(z)
    x = np.linspace(min(lnn),max(lnn),100)
    plt.plot(x,p(x),c='red')
    plt.errorbar(lnn,lntt,c='red',fmt='+')
    plt.show()

if __name__ == '__main__':
    main()
