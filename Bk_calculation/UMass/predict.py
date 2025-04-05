#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:53:42 2025

@author: kawall
"""
import numpy as np
import cmath
import matplotlib.pyplot as plt

#############################################################
## 
## this function takes as inputs
## x=radial position,y=vertical position in units of cm
## and outputs By form the transient normalized to the
## size of the transient at (x,y)=(0,0)
##
def transient(xin,yin):
    par0 = -1.0974293167
    par2 = +1.5724513306
    par4 = -1.0633523299
    par6 = +0.1774599020
    normy = 2.0*1.3011426132695032
    ntheta = 91
    theta_max = 2.0*np.pi*(31.5/360.0)
    dtheta = 2.0*theta_max/ntheta
    rplate = 3.748*3.50/2.902 #4.52033 ## 3.50*2.54/2.0   ## radius of kicker plates
    Bx = 0.0
    By = 0.0
    ## get contribution from first plate
    for i in range(0,ntheta):
        theta = -theta_max+(i+0.5)*dtheta
        arc = rplate*theta
        current=par0+par2*pow(arc,2.0)+par4*pow(arc,4.0)+par6*pow(arc,6.0)
        sx=rplate*np.cos(theta)
        sy=rplate*np.sin(theta)
        R = np.sqrt(pow(xin-sx,2.0)+pow(yin-sy,2.0))
        phi = cmath.phase(complex(xin-sx,yin-sy))
        Bx = Bx - current*np.sin(phi)/R
        By = By + current*np.cos(phi)/R
     
    ## get contribution from second plate by reversing sign of xin
    for i in range(0,ntheta):
        theta = -theta_max+(i+0.5)*dtheta
        arc = rplate*theta
        current=par0+par2*pow(arc,2.0)+par4*pow(arc,4.0)+par6*pow(arc,6.0)
        sx=rplate*np.cos(theta)
        sy=rplate*np.sin(theta)
        R = np.sqrt(pow(-xin-sx,2.0)+pow(yin-sy,2.0))
        phi = cmath.phase(complex(-xin-sx,yin-sy))
        Bx = Bx - current*np.sin(phi)/R
        By = By + current*np.cos(phi)/R
     
    return By/normy

##
##  end of function
############################################################

############################################################
##
##  try out the prediction at some points in (x,y)

## scan across x with y= 0
for j in range (0,33):
    x = -4.0+0.25*j
    y = 0.0
    BBy = transient(x,y)
    print("%5.3f %5.3f %8.5f" %(x,y,BBy))
   
   
## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=0.0 cm
sum = 0.0
for j in range (0,32):
    x = 0.00
    y = -1.6+(j+0.5)*0.1
    BBy = transient(x,y)
    print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=0.0",sum/32)

## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=1.75 cm
sum = 0.0
for j in range (0,32):
    x = 1.75
    y = -1.6+(j+0.5)*0.1
    BBy = transient(x,y)
    print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=1.75",sum/32)

## print full map with 0.5 mm resolution
#for x in np.arange(-45, 45.5, 0.5):
#    for y in np.arange(-45, 45.5, 0.5):
#        BBy = transient(0.1*x,0.1*y)
#        print("%5.1f %5.1f %8.5f" %(x,y,BBy))

print("B00: %5.3f" %(transient(0,0)))
print("B10: %5.3f" %(transient(1.75,0)))