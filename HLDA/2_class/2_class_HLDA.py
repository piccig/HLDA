#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import itertools
from scipy import linalg


N=6                   # Number of descriptors
Row_skip=0            # Skip <Row_skip> lines at the beginning (e.g. equilibration time)
Npoints=10            # Check convergence of the coefficintes taking <Npoints> chunks of the trajectory

fout=open('eigenvalues.dat','w')
fout2=open('eigenvectors.dat','w')

for j in range(1,Npoints+1):
    #state 1
    f1=open('1/COLVAR','r')
    data1= np.loadtxt(f1,usecols=range(1,N+1),skiprows=Row_skip)
    mu1=np.zeros(N)
    Nstep1 = int(len(data1)/Npoints)
    for i in range(0,N):
        mu1[i] = np.mean(data1[0:j*Nstep1,i])
    Cov1 = np.cov(np.transpose(data1))

    #state 2
    f2=open('2/COLVAR','r')
    data2= np.loadtxt(f2,usecols=range(1,N+1),skiprows=Row_skip)
    mu2=np.zeros(N)
    Nstep2 = int(len(data2)/Npoints)
    for i in range(0,N):
        mu2[i] = np.mean(data2[0:j*Nstep2,i])
    Cov2 = np.cov(np.transpose(data2))

#   CALCULATE POOLED COVARIANCE MATRIX AS SUM OF THE INVERSES OF SINGLE COVARINACES 
    Cov_inv = np.linalg.inv(Cov1) + np.linalg.inv(Cov2)

#   CALCULATE CENTROID OF THE DISTRIBUTIONS
    mutot=np.zeros(N)
    mutot = (mu1 + mu2)/2.0

#   DEFINE BETWEEN AND INVERSE WITHIN CLASS MATRICES
    between_class    = np.zeros((N,N))
    within_class_inv = np.zeros((N,N))

    between_class    = np.outer((mu1-mutot),(mu1-mutot)) + np.outer((mu2-mutot),(mu2-mutot))
    within_class_inv = Cov_inv

#   DEFINE MATRIX FOR EIGENDECOMPOSITION (S_W^-1 * S_B)
    tot_class = np.dot(within_class_inv,between_class)

#   EIGENDECOMPOSITION
    wt, vt = linalg.eig(tot_class)
    sidx = wt.argsort()[::-1] # sort according to the magnitude of the eigenvalues
    wt = wt[sidx]
    vt = vt[:,sidx]

    wt = np.real(wt)
    vt = np.real(vt)

#   THIS IS A PRINTING TRICK!
    vt = np.transpose(vt)

    wt = wt.tolist()
    vt = vt.tolist()

    time = j*Nstep1

    s = "  ".join(map(str, wt))
    txt=str(time)+'\t'+s+'\n'
    fout.write(txt)

    s = "  ".join(map(str, vt[0]))
    txt=str(time)+'\t'+s+'\n'
    fout2.write(txt)

print('Eigenvalues and eigenvector of Sigma_W^-1 * Sigma_B   ===> ', i)
for i in range(0,N):
    #print wt[i], vt[:,i]
    #print "%.3f,"*len(wt[i]) % tuple(wt[i]), "%.3f,"*len(vt[i]) % tuple(vt[i])
    print(wt[i], "%.3f,"*len(vt[i]) % tuple(vt[i]))

fout.close()
fout2.close()

