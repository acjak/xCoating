#!/usr/bin/env python
# encoding: utf-8
"""
cathode.py

Created by Anders Clemen Jakobsen on 2011-06-13.
Copyright (c) 2011 DTU Space, National Space Institute. All rights reserved.
"""

import matplotlib
#matplotlib.use('GTK3Agg')
import math,numpy,pylab
#import array 
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

class Plot:
    def __init__(self,glassLength,hDistance,mandrelRadius,boreDiameter,radiiOfLayers,d):
        pass
        
        self.glassLength = glassLength
        self.hDistance = hDistance
        self.mandrelRadius = mandrelRadius
        self.boreDiameter = boreDiameter
        self.radiiOfLayers = radiiOfLayers
        self.d = d
        
        self.drawOpticFromSide()
        self.drawOpticFromFront()
        
    def drawOpticFromSide(self):

        ax = pylab.figure().add_subplot(111)
        
        
        offset = 25.0 #mm
            
        parabStart = offset
        parabStop = offset + self.glassLength
        hyperbStart = offset + self.glassLength + self.hDistance
        hyperbStop = offset + 2*self.glassLength + self.hDistance
        
        bore = Rectangle(xy=[0.0, self.d-self.boreDiameter/2], width=20, height=self.boreDiameter, fill=False)    # empty rectangle
        mandrel = Rectangle(xy=[offset - 5.0, 20.0], width=hyperbStop-offset+5.0, height=25.0, fill=True)
        
        for i in range(len(self.radiiOfLayers)):
            parab = Line2D([parabStart,parabStop],[self.radiiOfLayers[i][0],self.radiiOfLayers[i][1]]) 
            hyperb = Line2D([hyperbStart,hyperbStop],[self.radiiOfLayers[i][3],self.radiiOfLayers[i][4]])                                   
            ax.add_line(parab)  
            ax.add_line(hyperb)      
            parab.set_color('blue')
            hyperb.set_color('red')
            parab.set_lw(1.0)
            hyperb.set_lw(1.0)
            if i == 0:
                 parab.set_color('black')
                 hyperb.set_color('black')
                 parab.set_lw(2.0)
                 hyperb.set_lw(2.0)
            
        bore.set_lw(2.0)
        bore.set_color('green')    
        mandrel.set_lw(1.0)
        mandrel.set_color('gray')
        ax.text(200,130,r'alpha1 $= 0.573$ degrees')
        ax.add_artist(mandrel)
        ax.add_artist(bore)
        pylab.xlim([0,500])
        pylab.ylim([0,150])
        pylab.savefig('OpticFromSide.ps',dpi=300,figsize=(6,6))
        #pylab.show()
        
    def drawOpticFromFront(self):
        ax2 = pylab.figure().add_subplot(111)
        bore = pylab.Circle((0,self.d),radius=self.boreDiameter/2,fill=False)
        ax2.add_artist(bore)
        #print len(self.radiiOfLayers)
        cover1 = Rectangle(xy=[0., 0.], width=-100.0, height=500.0, fill=True)
        cover2 = Rectangle(xy=[0., 0.], width=100.0, height=500.0, fill=True)
        t = matplotlib.transforms.Affine2D().rotate_deg(17.5)
        t_start = ax2.transData
        cover1.set_transform(t+t_start)
        t2 = matplotlib.transforms.Affine2D().rotate_deg(-17.5)
        cover2.set_transform(t2+t_start)


        for i in range(len(self.radiiOfLayers)):
            parab = pylab.Circle((0,0),radius=self.radiiOfLayers[i][0],fill=False,color='blue')
            hyperb = pylab.Circle((0,0),radius=self.radiiOfLayers[i][1],fill=False,color='red')
                               
            ax2.text(0,self.radiiOfLayers[i][0]-2.2,str(i),fontsize=9)
            #print self.radiiOfLayers[i][0]
            ax2.add_artist(parab)  
            ax2.add_artist(hyperb)
               


            #parab.set_lw(1.0)
            #hyperb.set_lw(1.0)
        
      
        

        cover1.set_lw(0.0)
        cover1.set_color('white')
        ax2.add_artist(cover1)
        cover2.set_lw(0.0)
        cover2.set_color('white')
        ax2.add_artist(cover2)
        ax2.add_artist(bore)

        parab = pylab.Circle((0,0),radius=self.radiiOfLayers[0][0],fill=True,color='gray')
        hyperb = pylab.Circle((0,0),radius=self.radiiOfLayers[0][0],fill=True,color='gray')
                                             
        ax2.add_artist(parab)  
        ax2.add_artist(hyperb)
        
        pylab.xlim([-60,60])
        pylab.ylim([20,120])
        pylab.savefig('OpticFromFront.ps',dpi=300,figsize=(6,6))
        pylab.show()
        
if __name__ == '__main__':
    pass