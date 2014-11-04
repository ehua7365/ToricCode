"""
torus.py
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import mwmatching

def main():
    test2()
##    for i in xrange(100):
##        for j in xrange(100):
##            (a1,b1) = test(i,j,10)
##            (a2,b2) = test(j,i,10)
##            if not (np.array_equal(a1,a2) or np.array_equal(b1,b2)):
##                print 'Error'
##                printLattice([],[],a1,b1,10)
##                print 'backwards'
##                printLattice([],[],a1,b1,10)
##                print (i,j)
##    test(0,0,10)
    generateErrors(10,0.1)
    for x in xrange(0,10):
        for y in xrange(0,100):
            if (x,y) != coordinates(index((x,y),10),10):
                print i
    

def plotCurves(Lmax,pmax,pstep,N):
    for L in xrange(2,Lmax):
        (x,y) = sampleCurve(L,pmax,pstep,N)
        z = np.polyfit(x,y,5)
        p = np.poly1d(z)
        xp = np.linspace(0,pmax,100)
        plt.plot(x,y,'.',xp,p(xp), '-')
        plt.ylim(0.0,1.0)
        plt.xlabel('Bit flip probability $p$')
        plt.ylabel('Information recovery probability $p_{recover}$')
        plt.title('Plot of $p_{recover}$ vs. $p$')
        plt.show()
        

def sampleCurve(L,pmax,pstep,N):
    pflip = np.empty(0,dtype=float)
    precover = np.empty(0,dtype=float)
    for p in xrange(pstep,pmax+pstep,pstep):
        n = 0
        for i in xrange(N):
            if (generateErrors(L,p)):
                n += 1
        np.append(pflip,p)
        np.append(precover,float(n)/N)
    return (pflip,precover)

def generateErrors(L,p):
    # Generate errors on each edge independently with probability p
    edgesX = np.less(np.random.rand(L,L),p) # Errors on horizontal edges
    edgesY = np.less(np.random.rand(L,L),p) # Errors on vertical edges
    (A,B) = findSyndromes(edgesX,edgesY,L)
    print 'lattice'
    printLattice(A,B,edgesX,edgesY,L)
    (pairsA,pairsB) = findPairs(A,B,edgesX,edgesY,L)
    correctErrors(edgesX,edgesY,pairsA,pairsB,L)
    (A,B) = findSyndromes(edgesX,edgesY,L)
    print 'correctedLattice'
    printLattice(A,B,edgesX,edgesY,L)
    return logicalX(edgesX,edgesY,L)&logicalZ(edgesX,edgesY,L)

def findPairs(A,B,edgesX,edgesY,L):
    # Generate the graphs for input into mwmatching algorithm
    nA,nB = len(A),len(B)
    graphEdgesA = [] # List of graph edges corresponding to path lengths A
    for i in xrange(nA-1):
        for j in xrange(i+1,nA):
            graphEdgesA.append((i,j,2*L-minDistance(A[i],A[j],L)))
    graphEdgesB = [] # List of graph edges corresponding to path lengths B
    for i in xrange(nB-1):
        for j in xrange(i+1,nB):
            graphEdgesB.append((i,j,2*L-minDistance(B[i],B[j],L)))

    # Feed graphs into mwmatching algorithm
    matchesA = mwmatching.maxWeightMatching(graphEdgesA)
    pairsA = []
    for i in xrange(len(matchesA)):
        p = [A[i],A[matchesA[i]]]
        p.sort()
        pairsA.append(p)
    # Remove duplicates
    pairsA = dict((x[0], x) for x in pairsA).values()
    
    matchesB = mwmatching.maxWeightMatching(graphEdgesB)
    pairsB = []
    for i in xrange(len(matchesB)):
        p = [B[i],B[matchesB[i]]]
        p.sort()
        pairsB.append(p)
    # Remove duplicates
    pairsB = dict((x[0], x) for x in pairsB).values()
    return (pairsA,pairsB)

def findSyndromes(edgesX,edgesY,L):
    A = np.empty(0,dtype=int) # Syndromes on vertices
    B = np.empty(0,dtype=int) # Syndromes on plaquettes
    # Find syndromes on vertices and plaquettes
    for x in xrange(L):
        for y in xrange(L):
            if (edgesX[x,y]^edgesY[x,y]
                ^edgesX[(x-1)%L,y]
                ^edgesY[x,(y-1)%L]):
                A = np.append(A,index((x,y),L))
            if (edgesX[x,y]^edgesY[x,y]
                ^edgesX[x,(y+1)%L]
                ^edgesY[(x+1)%L,y]):
                B = np.append(B,index((x,y),L))
    return (A,B)

def correctErrors(edgesX,edgesY,pairsA,pairsB,L):
    for pair in pairsA:
        drawShortestPathA(edgesX,edgesY,pair,L)
##    for pair in pairsB:
##        drawShortestPathB(edgesX,edgesY,pair,L)

def test2():
    L = 10
    x = np.zeros((L,L),dtype=bool)
    y = np.zeros((L,L),dtype=bool)
    x[2,5] = True
    (A,B) = findSyndromes(x,y,L)
    printLattice(A,B,x,y,L)
    print ''
    (pa,pb) = findPairs(A,B,x,y,L)
    print pa,pb
##    (A,B) = correctErrors(x,y,pa,pb,L)
    printLattice(A,B,x,y,L)
    

def test(i,j,L):
    edgesX = np.zeros((L,L),dtype=bool)
    edgesY = np.zeros((L,L),dtype=bool)
##    edgesX = np.matrix(
##        '0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0',dtype=bool)
##    edgesY = np.matrix(
##        '0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;\
##        0 0 0 0 0 0 0 0 0 0;',dtype=bool)
    drawShortestPathA(edgesX,edgesY,[i,j],L)
    return (edgesX,edgesY)

def drawShortestPathA(edgesX,edgesY,pair,L):
    (x0,y0) = coordinates(pair[0],L)
    (x1,y1) = coordinates(pair[1],L)
    dx = x1 - x0
    dy = y1 - y0
    if dx != 0:
        if dx < 0:
            temp = x0
            x0 = x1
            x1 = temp
        if abs(dx) <= L/2:
            xline = range(x0,x1)
        else:
            xline = range(x1,L)
            if x0 != 0:
                xline += range(0,x0)
##        print 'xline'
##        print xline
        for x in xline:
            edgesX[x,y0] ^= True # flip spins along x line
    if dy != 0:
        if dy < 0:
            temp = y0
            y0 = y1
            y1 = temp
        if abs(dy) <= L/2:
            yline = range(y0,y1)
        else:
            yline = range(y1,L)
            if y0 != 0:
                yline += (range(0,y0))
        for y in yline:
            edgesY[x0,y] ^= True # flip spins along y line
    
def drawShortestPathB(edgesX,edgesY,pair,L):
    pass

def logicalX(edgesX,edgesY,L):
    return True

def logicalZ(edgesX,edgesY,L):
    return True
    
def getMatchingPairs(adjacency):
    return True

def coordinates(i,L):
    return (i%L,i/L)

def index(coord,L):
    return coord[0]+L*coord[1]

def minDistanceAxis(xi,xj,L):
    x0 = abs(xi-xj)
    x1 = L-max(xi,xj)+min(xi,xj)
    return min(x0,x1)

def minDistance(i,j,L):
    (xi,yi) = coordinates(i,L)
    (xj,yj) = coordinates(j,L)
    return minDistanceAxis(xi,xj,L)+minDistanceAxis(yi,yj,L)
##    x0 = abs(xi-xj)
##    y0 = abs(yi-yj)
##    x1 = L-max(xi,xj)+min(xi,xj)
##    y1 = L-max(yi,yj)+min(yi,yj)
##    return min(x0,x1)+min(y0,y1)

def getEdgesTouchingVertex(v):
    n = len(v)/2
    L = int(math.sqrt(n))
    return np.array([v,n+v,(v/L+1)%L*L+v%L,(v%L+1)%L+v/L*L+n])

def getEdgesTouchingPlaquette(p):
    n = len(v)/2
    L = int(math.sqrt(n))
    return np.array([v,n+v,(v/L-1)%L*L+v%L,(v%L-1)%L+v/L*L+n])

def printLattice(A,B,edgesX,edgesY,L):
    for y in xrange(L):
        row1 = ''
        row2 = ''
        for x in xrange(L):
            if index((x,y),L) in A:
                row1 += 'X'
            else:
                row1 += ' '
            if edgesX[x,y]:
                row1 += '1'
            else:
                row1 += '0'
            if edgesY[x,y]:
                row2 += '1'
            else:
                row2 += '0'

            if index((x,y),L) in B:
                row2 += 'Z'
            else:
                row2 += ' '
        print row1
        print row2

def show(M):
    (n,m) = np.shape(M)
    for row in xrange(n):
        s = ''
        for col in xrange(m):
            if M[row,col]:
                s += '1'
            else:
                s += '0'
        print(s)

if __name__ == '__main__':
    main()
