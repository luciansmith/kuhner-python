# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:03:47 2016

@author: lpsmith
"""

from __future__ import division
import numpy as np
from scipy import optimize
from scipy.special import gamma
import matplotlib.pyplot as plt
from math import sqrt
from math import pi

#def model(p,x):
#    A,x1,sig1,B,x2,sig2 = p
#    return A*np.exp(-(x-x1)**2/sig1**2) + B*np.exp(-(x-x2)**2/sig2**2)

#def res(p,x,y):
#    return model(p,x) - y

#x = np.array([3963.67285156,  3964.49560547,  3965.31835938,  3966.14111328,  3966.96362305,
#         3967.78637695,  3968.60913086,  3969.43188477,  3970.25463867,  3971.07714844,
#         3971.89990234,  3972.72265625,  3973.54541016,  3974.36791992,  3975.19067383])
#y = np.array([1.75001533e-16,   2.15520995e-16,   2.85030769e-16,   4.10072843e-16, 7.17558032e-16,
#         1.27759917e-15,   1.57074192e-15,   1.40802933e-15, 1.45038722e-15,  1.55195653e-15,
#         1.09280316e-15,   4.96611341e-16, 2.68777266e-16,  1.87075114e-16,   1.64335999e-16])
#p0 = [1e-15,3968,2,1e-15,3972,2]
#p1,conv = optimize.leastsq(res,p0[:],args=(x,y))

#plt.plot(x,y,'+') # data
#fitted function
#plt.plot(np.arange(3962,3976,0.1),model(p1,np.arange(3962,3976,0.1)),'-')

def gaussian(pvec, val):
    mean, sigma, amp = pvec
    return abs(amp)*np.exp(-(val-mean)**2/sigma**2)
    
def tDistribution(pvec, val):
    mean, scale, amp, df = pvec
    df = abs(df)
    amp = abs(amp)
    return (gamma((df+1)/2)/sqrt(df*pi)) * (1 + (scale*((val-mean))**2)/df)**(-(df+1)/2)
    
def combined_gaussians(pvecs, val):
    ret = 0
    for index in range(0, int(len(pvecs)/3)):
        first = index*3
        pvec = pvecs[first], pvecs[first+1], pvecs[first+2]
        ret += gaussian(pvec, val)
    return ret
    
def combined_tdists(pvecs, val):
    ret = 0
    for index in range(0, int(len(pvecs)/4)):
        first = index*4
        pvec = pvecs[first], pvecs[first+1], pvecs[first+2], pvecs[first+3]
        ret += tDistribution(pvec, val)
    return ret
    
def getDataFrom(filename):
    datafile = open(filename, "r")
    xs = []
    ys = []
    for line in datafile:
        if (line.find("log2r") != -1):
            continue
        (x, y) = line.rstrip().split()
        if (str.isalpha(x[0])):
            continue
        xs.append(float(x))
        ys.append(float(y))
        #np.append(xs, x)
        #np.append(ys, y)
    return (np.array(xs), np.array(ys))
        
    
def resGauss(pvec,x,y):
    return combined_gaussians(pvec,x) - y

def resTdist(pvec,x,y):
    return combined_tdists(pvec,x) - y

def optimizeFile(filename, pvecs):
    xs, ys = getDataFrom(filename)
    optpvecs, conv = optimize.leastsq(resGauss, pvecs[:], args=(xs, ys))
    plt.plot(xs, ys, '.')
    plt.plot(np.arange(min(xs), max(xs),0.001), combined_gaussians(optpvecs, np.arange(min(xs), max(xs),0.001)), '-')
    for index in range(0, int(len(pvecs)/3)):
        first = index*3
        pvec = optpvecs[first], optpvecs[first+1], optpvecs[first+2]
        plt.plot(np.arange(min(xs), max(xs), 0.001), gaussian(pvec, np.arange(min(xs), max(xs), 0.001)), '--')
    plt.show()
    plt.close()
    for res in optpvecs:
        print res
    print conv

def optimizeFileTDist(filename, pvecs):
    xs, ys = getDataFrom(filename)
    optpvecs, conv = optimize.leastsq(resTdist, pvecs[:], args=(xs, ys))
    plt.plot(xs, ys, '.')
    plt.plot(np.arange(min(xs), max(xs),0.001), combined_tdists(optpvecs, np.arange(min(xs), max(xs),0.001)), '-')
    for index in range(0, int(len(optpvecs)/4)):
        first = index*4
        pvec = optpvecs[first], optpvecs[first+1], optpvecs[first+2], optpvecs[first+3]
        plt.plot(np.arange(min(xs), max(xs), 0.001), tDistribution(pvec, np.arange(min(xs), max(xs), 0.001)), '--')
    plt.show()
    plt.close()
    for res in optpvecs:
        print res
    print conv


optimizeFile("CN_rejoined_histograms/balanced_gain_hist_21-100000000.txt", [0.165, 0.0310917156, 4.18347237, 0.256, 0.10709327, 2.58297282, 0.354, 0.09286382, 0.42222656, 0.074, 0.09286382, 0.42222656])

optimizeFileTDist("CN_rejoined_histograms/balanced_gain_hist_21-100000000.txt", [0.171945606567, 117.683562259, 0.58297282, 0.0487337040807, 0.27068045802, 63.3855844071, 0.42222656, 0.116502162335, 0.354, 10.09286382, 0.42222656, 0.5, 0.074, 10.09286382, 0.42222656, 0.5])


    
#optimizeFile("CN-BAF_smoothed_histograms/all_hist_1.txt", [-0.63889836, 0.10917156, 0.18347237, -0.11649109, 0.10709327, 0.58297282, 0.19516854, 0.09286382, 3.42222656, 0.05252256, -0.67002912,  0.2604654])

#optimizeFile("ASCAT_smoothed_histograms/gain_hist__0-9000000000_all.txt", [-0.63889836, 0.10917156, 0.18347237, -0.11649109, 0.10709327, 0.58297282, 0.19516854, 0.09286382, 3.42222656])

#optimizeFile("ASCAT_smoothed_histograms/loss_hist__0-9000000000_all.txt", [-0.63889836, 0.10917156, 0.18347237, -0.11649109, 0.10709327, 0.58297282, -0.19516854, 0.09286382, 3.42222656])

#optimizeFile("ASCAT_smoothed_histograms/double_loss_hist__0-9000000000_.txt", [-0.63889836, 0.10917156, 0.18347237, -2.51649109, 0.10709327, 0.58297282, -1.9516854, 0.09286382, 3.42222656])

#optimizeFileTDist("ASCAT_smoothed_histograms/double_loss_hist__0-9000000000_.txt", [-0.63889836, 7.0917156, 0.18347237, 0.5, -2.21649109, 1.0709327, 0.58297282, 1, -1.9516854, 9.286382, 3.42222656, 1])

#optimizeFile("CN-BAF_smoothed_histograms/loss_hist_only__10-5000000_all.txt", [-0.106587123309, 0.0677821047683, 1.06303079888, -0.591018744029, -0.0984071240455, -0.493145131069, -0.255671670791, 0.126476044685, -1.009150511, -0.422971969928, 0.0866280356471, 0.95765122423, -0.596953397365, -0.277586365065, 0.835271786355])

#optimizeFile("CN-BAF_smoothed_histograms/gain_hist_only__10-5000000_all.txt", [0.106587123309, 0.0677821047683, 1.06303079888, 0.591018744029, -0.0984071240455, -0.493145131069, 0.255671670791, 0.126476044685, -1.009150511, 0.422971969928, 0.0866280356471, 0.95765122423, 0.596953397365, -0.277586365065, 0.835271786355])

#optimizeFile("CN-BAF_smoothed_histograms/double_loss_hist_only__10-5000000_all.txt", [-1.3889836, 0.10917156, 0.18347237, -2.61649109, 0.10709327, 0.58297282])

#optimizeFile("CN-BAF_smoothed_histograms/double_loss_hist_only__10-5000000_all.txt", [-1.3889836, 0.10917156, 0.18347237, -2.61649109, 0.10709327, 0.58297282, -2.61649109, 0.10709327, 0.58297282])

#optimizeFileTDist("CN-BAF_smoothed_histograms/double_loss_hist_only__10-5000000_all.txt", [-1.34842227851, 4.30309531052, 0.511360349936, 1, -2.49102730074, 7.008563427, 0.502895452949, 1])

#optimizeFileTDist("CN-BAF_smoothed_histograms/double_loss_hist_only__10-5000000_all.txt", [-1.34842227851, 4.30309531052, 0.511360349936, 10, -2.49102730074, 7.008563427, 0.502895452949, 10, -1.969102730074, 7.008563427, 0.502895452949, 10])

#plt.plot(np.arange(-2,2, 0.001), tDistribution([0, 1, 1, 1], np.arange(-2,2,0.001)), '--')
