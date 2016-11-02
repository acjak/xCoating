#!/usr/bin/env python
# encoding: utf-8

import sys
import os
import unittest
import matplotlib
#matplotlib.use('GTK3Agg')
import math,numpy,pylab
#import array 
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D


class Geometry:
    def __init__(self):
        #self.glassLength = glassLength
        #self.hDistance = hDistance
        #self.focalLength = focalLength
        #self.minimumRadius = minimumRadius
        #self.maximumRadius = maximumRadius
        #self.mandrelRadius = mandrelRadius
        #self.heightOfFirstLayer = 3.0 #mm
        #self.boreDiameter = boreDiameter
        #self.d = d
        #self.firstLayerOpening = firstLayerOpening
        
        self.radiiOfLayers = []
        self.boreOverlap = []
        self.areaOfLayer = []
        self.crossSectionArea = []
        self.effectiveGlass = []
        
        #self.buildOptic()
        #self.adjustMandrel()
        #self.calculateLayerOpeningToBore()
        #self.calculateAreaOfLayer()
        #self.calculateCrossSectionalAreaOfFullOptic()
        #self.effectiveGlassArea()
        
        #print len(self.crossSectionArea), len(self.boreOverlap)
        
        #for i in range(len(self.radiiOfLayers)):
            #print "Area of bore for layer ",i,": ",self.boreOverlap[i],". Cross sectional area for layer: ", self.crossSectionArea[i]
         #   print "Effective area of layer ",i,": ",self.effectiveGlass[i]," mm2"
        #print self.areaOfLayer
        
        
        #self.drawOpticFromSide()
        #self.drawOpticFromFront()
        
    def getGeometry(self,glassLength,hDistance,focalLength,minimumRadius,maximumRadius,mandrelRadius,boreDiameter,d,firstLayerOpening,glassThickness):
        self.glassLength = glassLength
        self.hDistance = hDistance
        self.focalLength = focalLength
        self.minimumRadius = minimumRadius
        self.maximumRadius = maximumRadius
        self.mandrelRadius = mandrelRadius
        #self.heightOfFirstLayer = 3.0 #mm
        self.boreDiameter = boreDiameter
        self.d = d
        self.firstLayerOpening = firstLayerOpening
        self.glassThickness = glassThickness
        
        self.radiiOfLayers.append(self.addMandrel()) #Add mandrel as 0th layer
        
        self.buildOptic()
        self.adjustMandrel()
        #self.calculateLayerOpeningToBore()
        #self.calculateAreaOfLayer()
        #self.calculateCrossSectionalAreaOfFullOptic()
        #self.effectiveGlassArea()
        
        return self.radiiOfLayers
        
    def getEffArea(self,d,radiiOfLayers,boreDiameter):
        
        self.d = d
        self.radiiOfLayers = radiiOfLayers
        self.boreDiameter = boreDiameter
        
        self.calculateLayerOpeningToBore()
        #print self.boreOverlap
        self.calculateAreaOfLayer()
        self.calculateCrossSectionalAreaOfFullOptic()
        self.effectiveGlassArea()
        
        return self.effectiveGlass
        
    def addMandrel(self):
        return [self.mandrelRadius,self.mandrelRadius,self.mandrelRadius,self.mandrelRadius,self.mandrelRadius]
        
    def adjustMandrel(self):
        R1,R2,R3,R4,R5 = self.radiiOfLayers[1]
        opening = self.firstLayerOpening
        self.radiiOfLayers[0] = [self.mandrelRadius,R2-opening,R3-opening,R4-opening,R5-opening]
        
    def buildOptic(self):
        self.i = 0
        while self.radiiOfLayers[self.i][0] <= self.d+self.boreDiameter/2:
            self.addLayer()
            self.i += 1
        
    def addLayer(self):
        R3 = self.radiiOfLayers[self.i][0]+self.glassThickness #R1 of layer underneath
        alpha = self.getAlpha(R3)
        #R2 = R3 + (self.hDistance/2)*math.tan(alpha)
        #R4 = R3 - (self.hDistance/2)*math.tan(3*alpha)
        R2 = R3 + (self.hDistance/2)*math.tan(alpha)
        R4 = R3 - (self.hDistance/2)*math.tan(3*alpha)
        R1 = R2 + math.tan(alpha)*self.glassLength
        R5 = R4 - math.tan(3*alpha)*self.glassLength
        self.radiiOfLayers.append([R1,R2,R3,R4,R5])
        
    def getAlpha(self,radius):
        alpha = math.atan(radius/self.focalLength)/4   
        return alpha
        
        
    def calculateAreaOfLayer(self):
        for i in range(len(self.radiiOfLayers)):
            R1,R2,R3,R4,R5 = self.radiiOfLayers[i]
            self.areaOfLayer.append(math.pi*math.sqrt(1+math.tan(self.getAlpha(R3))**2)*(R1**2+R3**2))
            #self.areaOfLayer[i] *= 
            
    def calculateCrossSectionalAreaOfFullOptic(self):
        for i in range(len(self.radiiOfLayers)):
            r0 = self.radiiOfLayers[i-1][0]
            r1 = self.radiiOfLayers[i][0]
            self.crossSectionArea.append(math.pi*(r1**2 - r0**2))
        return self.crossSectionArea

    def calculateLayerOpeningToBore(self):
        d = self.d #mm
        r = self.boreDiameter/2
        
        
        for i in range(len(self.radiiOfLayers)):
            R0 = self.radiiOfLayers[i-1][0] #Radius of layer below (R1)
            R1 = self.radiiOfLayers[i][0] #Radius of layer (R1)
            #print "R0 \t",R0,"R1 \t",R1
            if R0+r > d and R1+r > d and r+d > R1:
                A1_1 = r**2*math.acos((d**2+r**2-R1**2)/(2*d*r)) + R1**2*math.acos((d**2+R1**2-r**2)/(2*d*R1))
                A1_2 = 0.5*math.sqrt((-d+r+R1)*(d+r-R1)*(d-r+R1)*(d+R1+r))
                A1 = A1_1 - A1_2
        
                A0_1 = r**2*math.acos((d**2+r**2-R0**2)/(2*d*r)) + R0**2*math.acos((d**2+R0**2-r**2)/(2*d*R0))
                A0_2 = 0.5*math.sqrt((-d+r+R0)*(d+r-R0)*(d-r+R0)*(d+R0+r))
                A0 = A0_1 - A0_2 
        
                A = A1 - A0
            
                self.boreOverlap.append(A)
                #print i, A, "mm2 overlap", "  Case 1"
            
            elif d-r < R1 and d-r > R0:
                A1_1 = r**2*math.acos((d**2+r**2-R1**2)/(2*d*r)) + R1**2*math.acos((d**2+R1**2-r**2)/(2*d*R1))
                A1_2 = 0.5*math.sqrt((-d+r+R1)*(d+r-R1)*(d-r+R1)*(d+R1+r))
                A = A1_1 - A1_2
            
                self.boreOverlap.append(A)
                
                #print i, A, "mm2 overlap", "  Case 2"
            
            elif d+r > R0 and d+r < R1:
                A_bore = math.pi*r**2
                A0_1 = r**2*math.acos((d**2+r**2-R0**2)/(2*d*r)) + R0**2*math.acos((d**2+R0**2-r**2)/(2*d*R0))
                A0_2 = 0.5*math.sqrt((-d+r+R0)*(d+r-R0)*(d-r+R0)*(d+R0+r))
                A0 = A0_1 - A0_2
            
                A = A_bore - A0
            
            
                self.boreOverlap.append(A)
                #print i, A, "mm2 overlap", "  Case 3"
            
            else:
                pass
                self.boreOverlap.append(0)
                #print i, "No overlap"

        return self.boreOverlap
   
        
    def effectiveGlassArea(self):
        for i in range(len(self.radiiOfLayers)):
            self.effectiveGlass.append(self.areaOfLayer[i]*self.boreOverlap[i]/self.crossSectionArea[i])
            
        
    
    
if __name__ == '__main__':
    Geometry()