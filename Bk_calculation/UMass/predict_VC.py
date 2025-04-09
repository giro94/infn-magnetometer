#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 16:47:54 2025

@author: kawall
"""
import numpy as np
import cmath
import matplotlib.pyplot as plt


############################################################
##
## read in the estimated current distributionn in VC
##
vc = np.loadtxt('vc_currents.txt',delimiter='\t',skiprows=1,usecols=(0,1,5))

vcx = vc[0:,0]   ## source x in umass coordinates
vcy = vc[0:,1]   ## source y in umass coordinates
vcj = vc[0:,2]   ## vacuum chamber current density

#############################################################
## 
## this function takes as inputs FNAL positions
## x=radial position,y=vertical position in units of cm
## and outputs By from the transient originating from
## the vacuum chamber+kicker cage, normalized to the
## observed size of vc_transient(x=0,y=0)=-35.0
##
def vc_transient(xin_raw,yin_raw,vcx_in,vcy_in,vcj_in):
    RAT  = 3.50/2.902   ## Ratio of FNAL/UMass dimensions
    xin = xin_raw/RAT
    yin = yin_raw/RAT
    nsteps = len(vcx_in)

    Bx = 0.0
    By = 0.0
    scale_vc = 1.3260

    normy=1.0
    for i in range(0,nsteps):
        current=vcj_in[i]
        sx=vcx_in[i]
        sy=vcy_in[i]
        current=vcj_in[i]

        ## get contribution from upper aluminum plate
        R = np.sqrt(pow(xin-sx,2.0)+pow(yin-sy,2.0))
        phi = cmath.phase(complex(xin-sx,yin-sy))
        Bx = Bx - current*np.sin(phi)/R
        By = By + current*np.cos(phi)/R

        ## get contribution from lower aluminum plate by reversing sign of yin
        sy=-vcy_in[i]
        R = np.sqrt(pow(xin-sx,2.0)+pow(yin-sy,2.0))
        phi = cmath.phase(complex(xin-sx,yin-sy))
        Bx = Bx + current*np.sin(phi)/R
        By = By + current*np.cos(phi)/R

    return -(By/normy)*scale_vc

##
##  end of function
############################################################
#############################################################
## 
## this function takes as inputs FNAL positions
## x=radial position, y=vertical position in units of cm
## and outputs By from the transient normalized to make the
## size of the kick+transient at plates_transient(x=0,y=0) = 189
## to match the data

def plates_transient(xin_raw,yin_raw):
    RAT  = 3.50/2.902   ## Ratio of FNAL/UMass dimensions
    xin = xin_raw/RAT
    yin = yin_raw/RAT

    par0 = -1.101319  ## use this value to predict transient, not kick
    ## par0 = +1.175  ## use this value to predict kick+transient
    par2 = +2.275000
    par4 = -2.24500
    par6 = +0.546500
    normy = (198.0/101.88681)
    scale_plates2 = 0.95449

    ntheta = 91
    theta_max = 2.0*np.pi*(31.5/360.0)
    dtheta = 2.0*theta_max/ntheta
    rplate = 3.748 #3.748   ## radius of UMass kicker plates
    Bx = 0.0
    By = 0.0
    ## get contribution from first plate
    for i in range(0,ntheta):
        theta = -theta_max+(i+0.5)*dtheta
        arc = rplate*theta
        current=par0+par2*pow(arc,2.0)+par4*pow(arc,4.0)+par6*pow(arc,6.0)
        sx=rplate*np.cos(theta)
        sy=rplate*np.sin(theta)

        ## get contribution from first plate
        R = np.sqrt(pow(xin-sx,2.0)+pow(yin-sy,2.0))
        phi = cmath.phase(complex(xin-sx,yin-sy))
        Bx = Bx - current*np.sin(phi)/R
        By = By + current*np.cos(phi)/R

        ## get contribution from second plate by reversing sign of xin
        R = np.sqrt(pow(-xin-sx,2.0)+pow(yin-sy,2.0))
        phi = cmath.phase(complex(-xin-sx,yin-sy))
        Bx = Bx + current*np.sin(phi)/R
        By = By + current*np.cos(phi)/R

    return -By*normy*scale_plates2

##
##  end of function
############################################################

def transient(x,y):
    return (plates_transient(x,y)+vc_transient(x,y,vcx,vcy,vcj))/(plates_transient(0,0)+vc_transient(0,0,vcx,vcy,vcj))






## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=0.0 cm
sum = 0.0
for j in range (0,321):
    x = 0.00
    y = -1.6+j*0.01
    BBy = transient(x,y)
    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=0.0",sum/321)

## scan over UMass crystal from y=-1.45 to y=+1.45 cm at x=2.9 cm
sum = 0.0
for j in range (0,291):
    x = 0.00
    y = -1.45+j*0.01
    BBy = transient(x,y)
    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("UMass crystal at x=0.0",sum/291)

## scan over UMass crystal from y=-1.45 to y=+1.45 cm at x=2.9 cm
#sum = 0.0
#for j in range (0,59):
#    x = 0.00
#    y = -1.45+j*(2.9/58.0)
#    BBy = transient(x,y)
#    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
#    sum=sum+BBy
#   
#print("UMass crystal at x=0.0 (59 steps)",sum/59)

## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=1.75 cm
sum = 0.0
for j in range (0,321):
    x = 1.75
    y = -1.6+j*0.01
    BBy = transient(x,y)
    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=1.75",sum/321)


## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=1.75 cm
sum = 0.0
for j in range (0,321):
    x = 0.2
    y = -1.6+j*0.01
    BBy = transient(x,y)
    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=0.2",sum/321)

## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=1.75 cm
sum = 0.0
for j in range (0,321):
    x = 0.0
    y = -1.6+j*0.01 +0.2
    BBy = transient(x,y)
    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=0.0, y+0.2",sum/321)

## scan over INFN crystal from y=-1.6 to y=+1.6 cm at x=1.75 cm
sum = 0.0
for j in range (0,321):
    x = 0.0
    y = -1.6+j*0.01 +0.5
    BBy = transient(x,y)
    #print("%5.3f %5.3f %8.5f" %(x,y,BBy))
    sum=sum+BBy
   
print("INFN crystal at x=0.0, y+0.5",sum/321)

print("B00: %5.3f" %(transient(0,0)))
print("B10: %5.3f" %(transient(1.75,0)))