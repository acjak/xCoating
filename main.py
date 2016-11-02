#!/usr/bin/env python
# encoding: utf-8


from geometry import Geometry
from plot import Plot
#from reflectivity import Reflectivity
#from coating import Coating
import pidly
import numpy,math
#from scipy import convolve,square
import pylab
from matplotlib.pylab import plot,show,xlabel,ylabel,legend,title,subplots_adjust,savefig,cla,hist
import matplotlib.pyplot as plt
from time import localtime, strftime
import os

class Main:
    def __init__(self):

        self.glassThickness = 0.21 #mm
        #CAST
        self.doCAST()

        #IAXO
        #self.doIAXO()

        self.N_recipes =4.0
        self.HPD = 48. #60 arc seconds

        #Coating parameters
        #self.matCombArray = [['Pt_llnl_cxro','B4C_llnl_cxro'],['W_llnl_cxro','B4C_llnl_cxro'],['W_llnl_cxro','Si_llnl_cxro'],['Ni.93V.07_llnl_cxro','B4C_llnl_cxro'],['Ni.93V.07_llnl_cxro','a-C_llnl_cxro'],['Pt_llnl_cxro','a-C_llnl_cxro']]
        #self.matCombArray = [['Pt_llnl_cxro','B4C_llnl_cxro'],['W_llnl_cxro','B4C_llnl_cxro'],['Ni.93V.07_llnl_cxro','B4C_llnl_cxro']]
        self.matCombArray = [['Pt_llnl_cxro','a-C_llnl_cxro']]

        self.subs = 'Si_llnl_cxro'
        #self.mat1 = 'Pt_llnl_cxro'
        #self.mat2 = 'B4C_llnl_cxro'
        self.subsRoughness = 3.5
        self.mat1Roughness = 5.0
        self.mat2Roughness = 5.0
        self.n_min = 1
        self.n_max = 10.1
        self.gamma_min = 0.3
        self.gamma_max = 0.7
        self.gradedD = True

        self.slope_min = 30.0
        self.slope_max = 400.0

        self.nStepSize = 1.0
        self.dStepSize = 5.0
        self.gammaStepSize = 0.05

        self.energyValues = 100.
        self.eMin = 0.1
        self.eMax = 10.0

        #Single layer:
        #self.SLMat_array = ['Pt_llnl_cxro','W_llnl_cxro','Ni.93V.07_llnl_cxro','Ir_llnl_cxro']
        self.SLMat_array = ['Ni.93V.07_llnl_cxro']
        self.D = 500
        self.matR = 5.0

        if self.gamma_min != self.gamma_max:
            self.gamma = numpy.arange(self.gamma_min,self.gamma_max,self.gammaStepSize)
        else:
            self.gamma = [self.gamma_min]

        if self.n_min != self.n_max:
            self.n = numpy.arange(self.n_min,self.n_max,self.nStepSize)
        else:
            self.n = [self.n_min]


        #self.slope_d_min = [self.slope_min]
        #self.slope_d_max = [self.slope_max]
        #self.slope_d_max = numpy.arange(self.slope_min,self.slope_max,self.dStepSize)


        if self.slope_min == self.slope_max:
            self.slope_d_min = [self.slope_min]
            self.slope_d_max = [self.slope_max]
        else:
            self.slope_d_min = numpy.arange(self.slope_min,self.slope_max,self.dStepSize)
            self.slope_d_max = numpy.arange(self.slope_min,self.slope_max,self.dStepSize)



        self.time = strftime("%d%m%y-%H%M", localtime())
        self.directory = 'output/'+self.time
        print self.directory





        self.geo = Geometry()
        self.getGeometry(self.focalLength)
        print self.radiiOfLayers
        #self.drawOptic(self.focalLength)
        self.getGeometricArea()

        #self.N_per_recipe = int(float(len(self.radiiOfLayers))/float(self.N_recipes))
        self.N_per_recipe = 4



        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        numpy.savetxt('output/%s/radiiOfLayers.txt' % (self.time),self.radiiOfLayers,delimiter='\t')
        with open("output/%s/anglesOfLayers.txt" % (self.time), "a") as myfile:
            myfile.write("Layer \t Area (mm2) \t Angle (degrees) \t Radius (mm)")

        self.ML_mats = []

        for i in range(len(self.matCombArray)):
            self.ML_mats.append(self.matCombArray[i][0][:-10] + '-' + self.matCombArray[i][1][:-10])

        for i in range(len(self.matCombArray)):
            os.makedirs(self.directory+'/'+str(self.ML_mats[i]))

        self.SL_mats = []
        for i in range(len(self.SLMat_array)):
            self.SL_mats.append(self.SLMat_array[i][:-10])
        for i in range(len(self.SL_mats)):
            os.makedirs(self.directory+'/'+str(self.SL_mats[i]))

        self.optimizeForCoating()

    def doCAST(self):
        self.glassLength = 225.0 #mm
        self.hDistance = 4.0 #mm
        self.focalLength = 1500.0 #mm
        self.minimumRadius = 53.0 #mm
        self.maximumRadius = 100.0 #mm
        self.mandrelRadius = 60.5 #mm
        self.boreDiameter = 43.0 #mm
        self.d = 83.0 #mm
        self.firstLayerOpening = 2.5 #mm

    def doIAXO(self):
        self.glassLength = 225.0 #mm
        self.hDistance = 2.0 #mm
        self.focalLength = 7000.0  # mm
        self.minimumRadius = 50.0 #mm
        self.maximumRadius = 350.0 #mm
        self.mandrelRadius = 53.0 #mm
        self.boreDiameter = 600.0 #mm
        self.d = 0.0 #mm
        self.firstLayerOpening = 2.5 #mm

    def optimizeForFocallength(self):
        ##############
        self.focalLength = [4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
        #self.focalLength = [5000.0,7000.0]

        #MLA_effect = numpy.zeros((len(self.focalLength),self.energyValues+1))
        ML_Int = numpy.zeros((len(self.focalLength)))
        SL_Int = numpy.zeros((len(self.focalLength)))

        ML_Int_w_spot = numpy.zeros((len(self.focalLength)))
        SL_Int_w_spot = numpy.zeros((len(self.focalLength)))

        intTotalThroughput = numpy.zeros((len(self.focalLength)))
        spotsizeArray = numpy.zeros((len(self.focalLength)))

        cla()
        fig = plt.figure()
        fig2 = plt.figure()
        fig3 = plt.figure()
        fig4 = plt.figure()
        fig5 = plt.figure()
        fig6 = plt.figure()
        fig7 = plt.figure()
        fig8 = plt.figure()
        fig9 = plt.figure()
        #ba = fig.add_subplot(111)
        #bc = fig3.add_subplot(111)
        #SL_TP_plot = fig2.add_subplot(111)
        #ML_TP_plot = fig4.add_subplot(111)
        #SL_TP_plot_spot = fig5.add_subplot(111)
        #ML_TP_plot_spot = fig6.add_subplot(111)

        plot1a = fig.add_subplot(111)
        plot1b = fig2.add_subplot(111)
        plot2a = fig3.add_subplot(111)
        plot2b = fig4.add_subplot(111)
        plot3 = fig5.add_subplot(111)
        plot4 = fig6.add_subplot(111)
        plot5 = fig7.add_subplot(111)
        plot6 = fig8.add_subplot(111)
        plot7 = fig9.add_subplot(111)

        plot7_2 = plot7.twinx()

        boreArea = math.pi*(0.5*self.boreDiameter/10)**2

        dE = (self.eMax-self.eMin)/(self.energyValues+1)
        energy = numpy.arange(self.eMin,self.eMax,dE)
        aspectrum = self.getAxionSpectrum(energy)

        qe = numpy.loadtxt('Eff_mM_No_Strongback.txt', delimiter = '\t')

        det_qe = numpy.interp(energy,qe[:,0],qe[:,1])

        for p in range(len(self.focalLength)):
            print "\n\n\n\nFocal length = ", self.focalLength[p]
            self.geo = Geometry()

            self.getGeometry(self.focalLength[p])

            #spotSize = self.getSpotSize50(self.focalLength[p])
            spotSize = self.getSpotSize50(self.focalLength[p])
            spotsizeArray[p] = math.sqrt(spotSize)

            self.N_per_recipe = int(float(len(self.radiiOfLayers))/float(self.N_recipes))

            self.getGeometricArea()
            #numpy.savetxt('output/%s/radiiOfLayers_%s.txt' % (self.time,self.focalLength[p]),self.radiiOfLayers,delimiter='\t')

            print "Opening area received"

            #MLA_effect,ML_Int[p] = self.calculateReflectivityInIDL(self.focalLength[p])
            #SLAeffective,SL_Int[p] = self.calculateSLReflectivityIDL(self.focalLength[p])

            #MLA_throughput = [(100*i)/boreArea for i in MLA_effect[0]]
            #SLA_throughput = [(100*i)/boreArea for i in SLAeffective[0]]

            #MLA_throughput_w_spotsize = [i/math.sqrt(spotSize) for i in MLA_throughput]
            #SLA_throughput_w_spotsize = [i/math.sqrt(spotSize) for i in SLA_throughput]

            #ML_tofile = []
            #SL_tofile = []
            #for k in range(len(self.idl.energy)):
            #    ML_tofile.append([self.idl.energy[k],MLA_effect[0][k]])
            #    SL_tofile.append([self.idl.energy[k],SLAeffective[0][k]])

            #numpy.savetxt('output/%s/ML_effective_area_%s.txt' % (self.time,self.focalLength[p]),ML_tofile,delimiter='\t')
            #numpy.savetxt('output/%s/SL_effective_area_%s.txt' % (self.time,self.focalLength[p]),SL_tofile,delimiter='\t')

            #ML_Int_w_spot[p] = (ML_Int[p]/boreArea)/math.sqrt(spotSize)
            #SL_Int_w_spot[p] = (SL_Int[p]/boreArea)/math.sqrt(spotSize)

            #print len(MLA_effect),len(self.idl.energy)
            #print MLA_effect[0]

            refl, refl_w_detector, eff_area, eff_area_w_det, totalcrossarea = self.makePlotsForLOI(self.focalLength[p])

            throughput = [(100*i)/boreArea for i in refl] # Optics throughput
            # Optics eff. area
            throughput_det = [(100*i)/boreArea for i in refl_w_detector] # Optics + detector throughput
            # Optics + detector eff. area

            eff_area_w_det_secondary = [i*j for i,j in zip(det_qe,eff_area)]

            total_throughput = [i*j*totalcrossarea for i,j in zip(throughput_det,aspectrum)]
            total_throughput_total = [i for i in total_throughput]
            #throughput_det * aspectrum# Total axion throughput

            total_throughput_w_spot = [i/math.sqrt(spotSize) for i in total_throughput_total] # Total axion throughput / sqrt(spot)

            intTotalThroughput[p] = sum([n*dE for n in total_throughput_w_spot])

            print str(totalcrossarea) + " m2"

            ML_tofile = []
            for k in range(len(energy)):
                ML_tofile.append([energy[k],throughput[k]])
            numpy.savetxt('output/%s/Throughput_vs_energy_%s.txt' % (self.time,self.focalLength[p]),ML_tofile,delimiter='\t')

            plot1a.plot(energy,throughput,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot1b.plot(energy,eff_area,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot2a.plot(energy,throughput_det,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot2b.plot(energy,eff_area_w_det_secondary,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot3.plot(energy,total_throughput_total,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot4.plot(energy,total_throughput_w_spot,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot5.plot(energy,total_throughput_w_spot,label='Focal length: %s m' % (self.focalLength[p]/1000))

            plot7.plot(energy,throughput,label='Focal length: %s m' % (self.focalLength[p]/1000))
            plot7_2.plot(energy,eff_area,label='Focal length: %s m' % (self.focalLength[p]/1000))


            #ba.plot(self.idl.energy,MLA_effect[0],label='Focal length: %s' % (self.focalLength[p]))
            #bc.plot(self.idl.energy,SLAeffective[0],label='Focal length: %s' % (self.focalLength[p]))
            #ML_TP_plot.plot(self.idl.energy,MLA_throughput,label='Focal length: %s' % (self.focalLength[p]))
            #SL_TP_plot.plot(self.idl.energy,SLA_throughput,label='Focal length: %s' % (self.focalLength[p]))
            #ML_TP_plot_spot.plot(self.idl.energy,MLA_throughput_w_spotsize,label='Focal length: %s' % (self.focalLength[p]))
            #SL_TP_plot_spot.plot(self.idl.energy,SLA_throughput_w_spotsize,label='Focal length: %s' % (self.focalLength[p]))

            #ba.plot(self.idl.energy,MLA_effect_w_spotsize_m,label='Focal length: %s' % (self.focalLength[p]))
            #bc.plot(self.idl.energy,SLA_effect_w_spotsize_m,label='Focal length: %s' % (self.focalLength[p]))
            #self.drawOptic(self.focalLength[p])


        print "test"

        plot1a.set_xlabel('Energy [keV]')
        plot1a.set_ylabel('Throughput [%]')
        #ba.set_title('Effective reflectivity per area comparison')
        plot1a.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        fig.savefig('output/%s/Throughput.eps' % (self.time))
        cla()

        print "test"

        plot1b.set_xlabel('Energy [keV]')
        plot1b.set_ylabel('Effective area [cm2]')
        #ba.set_title('Effective reflectivity per area comparison')
        plot1b.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        fig2.savefig('output/%s/Effective_area.eps' % (self.time))
        cla()

        plot2a.set_xlabel('Energy [keV]')
        plot2a.set_ylabel('Throughput * detector efficiency [%]')
        #ba.set_title('Effective reflectivity per area comparison')
        plot2a.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        fig3.savefig('output/%s/Throughput_w_detector.eps' % (self.time))
        cla()

        plot2b.set_xlabel('Energy [keV]')
        plot2b.set_ylabel('Effective area [cm2] * detector efficiency')
        #ba.set_title('Effective reflectivity per area comparison')
        plot2b.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        fig4.savefig('output/%s/Effective_area_w_detector.eps' % (self.time))
        cla()

        plot3.set_xlabel('Energy [keV]',fontsize='15')
        #plot3.set_ylabel(r'$\frac{d\Phi}{dE} (10^{20} \textnormal{keV}^{-1} \textnormal{year}^{-1} \textnormal{m}^{-2})$')

        plot3.set_ylabel('DAF [Arb. units]',fontsize='15')
        #ba.set_title('Effective reflectivity per area comparison')
        #plot3.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        rightarrow = dict(boxstyle="rarrow,pad=0.3", fc="white", ec="black", lw=1)
        leftarrow = dict(boxstyle="larrow,pad=0.3", fc="white", ec="black", lw=1)
        t = plot3.text(3.5, 10, "4 m focal length", ha="center", va="center", rotation=45,
            size=12,
            bbox=rightarrow)
        t2 = plot3.text(6.25, 35, "10 m focal length", ha="center", va="center", rotation=25,
            size=12,
            bbox=leftarrow)
        fig5.savefig('output/%s/Total_throughput.eps' % (self.time))
        cla()

        plot4.set_xlabel('Energy [keV]',fontsize='15')
        plot4.set_ylabel(r'$\mathrm{DAF} * a^{-\frac{1}{2}} \left[\mathrm{Arb. units }\right]$',fontsize='15')
        #ba.set_title('Effective reflectivity per area comparison')
        #plot4.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        fig6.savefig('output/%s/Total_throughput_w_spot.eps' % (self.time))
        cla()

        plot5.set_xlabel('Energy [keV]')
        plot5.set_ylabel(r'$\frac{d\Phi}{dE}$')
        #ba.set_title('Effective reflectivity per area comparison')
        #fig7.savefig('output/%s/Total_throughput_w_spot2.png' % (self.time))
        plot5.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        cla()

        for i in range(len(self.focalLength)):
            self.focalLength[i] = self.focalLength[i]/1000

        plot6.plot(self.focalLength,intTotalThroughput,'bo',self.focalLength,intTotalThroughput,'b',ms=8)
        plot6.set_xlabel('Focal length [m]',fontsize='15')
        plot6.set_ylabel(r'$\mathrm{DAF} * a^{-\frac{1}{2}} \left[\mathrm{Arb. units }\right]$',color='b',fontsize='15')
        #ba.set_title('Effective reflectivity per area comparison')
        plot6.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        for tl in plot6.get_yticklabels():
            tl.set_color('b')
        y2 = plot6.twinx()
        y2.plot(self.focalLength,spotsizeArray,'rs',self.focalLength,spotsizeArray,'r--',ms=8)
        y2.set_ylabel(r'$a^{\frac{1}{2}} [ \mathrm{cm} ]$',color='r',fontsize='15')
        for tl in y2.get_yticklabels():
            tl.set_color('r')
        fig8.savefig('output/%s/int_total_throughput_w_spot.eps' % (self.time))

        cla()

        plot7.set_xlabel('Energy [keV]',fontsize='15')
        plot7.set_ylabel('Throughput [%]',fontsize='15')
        #plot7.set_ylim([0,100])
        plot7_2.set_ylabel('Effective area [cm2]',fontsize='15')
        #for tl in plot7.get_yticklabels():
        #    tl.set_color('b')
        #for tl in plot7_2.get_yticklabels():
        #    tl.set_color('r')


        rightarrow = dict(boxstyle="rarrow,pad=0.3", fc="white", ec="black", lw=1)
        leftarrow = dict(boxstyle="larrow,pad=0.3", fc="white", ec="black", lw=1)
        t = plot7.text(2.5, 30, "4 m focal length", ha="center", va="center", rotation=45,
            size=12,
            bbox=rightarrow)
        t2 = plot7.text(7.5, 75, "10 m focal length", ha="center", va="center", rotation=25,
            size=12,
            bbox=leftarrow)
        #ba.set_title('Effective reflectivity per area comparison')
        #plot1a.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        fig9.savefig('output/%s/Throughput_and_eff_area.eps' % (self.time))
        cla()



        print "test"

        # ba.set_xlabel('Energy [keV]')
        # ba.set_ylabel('Effective area [cm2]')
        # #ba.set_title('Effective reflectivity per area comparison')
        # ba.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        # fig.savefig('output/%s/Effective_area_ML.png' % (self.time))
        # cla()
        # bc.set_xlabel('Energy [keV]')
        # bc.set_ylabel('Effective area [cm2]')
        # #bc.set_title('Effective reflectivity per area comparison')
        # bc.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        # fig3.savefig('output/%s/Effective_area_SL.png' % (self.time))
        # cla()
        # SL_TP_plot.set_xlabel('Energy [keV]')
        # SL_TP_plot.set_ylabel('Throughput [%]')
        # #SL_TP_plot.set_title('Effective reflectivity per area comparison')
        # SL_TP_plot.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        # fig2.savefig('output/%s/Throughput_SL.png' % (self.time))
        # cla()
        # ML_TP_plot.set_xlabel('Energy [keV]')
        # ML_TP_plot.set_ylabel('Throughput [%]')
        # #ML_TP_plot.set_title('Effective reflectivity per area comparison')
        # ML_TP_plot.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        # fig4.savefig('output/%s/Throughput_ML.png' % (self.time))
        # cla()
        # SL_TP_plot_spot.set_xlabel('Energy [keV]')
        # SL_TP_plot_spot.set_ylabel('Throughput / sqrt(spot)')
        # #SL_TP_plot_spot.set_title('Effective reflectivity per area comparison')
        # SL_TP_plot_spot.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        # fig5.savefig('output/%s/Throughput_spot_SL.png' % (self.time))
        # cla()
        # ML_TP_plot_spot.set_xlabel('Energy [keV]')
        # ML_TP_plot_spot.set_ylabel('Throughput / sqrt(spot)')
        # #ML_TP_plot_spot.set_title('Effective reflectivity per area comparison')
        # ML_TP_plot_spot.legend(loc='upper center', bbox_to_anchor=(0.8,1.0))
        # fig6.savefig('output/%s/Throughput_spot_ML.png' % (self.time))
        # cla()


        #MLAeffective,MLIntRefl = self.calculateReflectivityInIDL()
        #SLAeffective,SLIntRefl = self.calculateSLReflectivityIDL()
        #self.removeEndOfString()

        #self.plotEffectiveArea(MLAeffective,SLAeffective)
    #self.plotIntReflComp(ML_Int_w_spot,SL_Int_w_spot,self.focalLength)
        #self.plotIntReflComp(MLAeffective,SLAeffective)
        #self.saveEffAreaToFile(MLAeffective,SLAeffective)


        #for i in range(len(self.radiiOfLayers)):
        #    print self.radiiOfLayers[i]
        #print len(self.radiiOfLayers)
        #self.drawOptic()


    def optimizeForCoating(self):
        self.idl = pidly.IDL()#long_delay=0.0000001)
        self.idl.pro('.r imd')

        #####################!!!!!!!!!!!!!!!!##################
        MLAeffective,MLIntRefl = self.calculateReflectivityInIDL(self.focalLength)
        self.removeEndOfString()

        self.saveEffAreaToFile(MLAeffective)
        #SLAeffective,SLIntRefl = self.calculateSLReflectivityIDL()
        self.plotEffectiveArea(MLAeffective)#,SLAeffective)
        #self.plotIntComp(self.getIntArray(MLIntRefl,SLIntRefl))
        #self.plotReflComp(self.getAreaArray(MLAeffective,SLAeffective))
        #self.saveEffAreaToFile(MLAeffective,SLAeffective)





    def onlyDoGeometry(self):
        pass

    def getSpotSize50(self,focalLength):
        EED = 2*(focalLength/10)*math.tan(math.radians(0.5*self.HPD/(60.*60.)))
        sun_EED = 2*(focalLength/10)*math.tan(math.radians(0.5*3.16/60.))

        EED_w_sun = math.sqrt(EED**2 + sun_EED**2)

        spotArea = math.pi*(0.5*EED_w_sun)**2

        #math.sqrt((math.pi*(focalLength*math.tan(math.radians(1./self.HPD)))**2)/10**6)

        return spotArea

    def getSpotSize80(self,focalLength):
        EED = 2*(focalLength/10)*math.tan(math.radians(0.5*self.HPD/(60.*60.)))
        sun_EED = 2*(focalLength/10)*math.tan(math.radians(0.5*3.61/60.))

        EED_w_sun = math.sqrt(EED**2 + sun_EED**2)

        spotArea = math.pi*(0.5*EED_w_sun)**2

        #math.sqrt((math.pi*(focalLength*math.tan(math.radians(1./self.HPD)))**2)/10**6)

        return spotArea

    def getGeometry(self,focalLength):
        self.radiiOfLayers = self.geo.getGeometry(self.glassLength,self.hDistance,focalLength,self.minimumRadius,self.maximumRadius,self.mandrelRadius,self.boreDiameter,self.d,self.firstLayerOpening,self.glassThickness)
        #print self.radiiOfLayers

    def getGeometricArea(self):
        self.geoArea = self.geo.getEffArea(self.d,self.radiiOfLayers,self.boreDiameter)
        self.crossArea = self.geo.calculateLayerOpeningToBore()
        #self.crossArea = self.geo.calculateCrossSectionalAreaOfFullOptic()
        #print self.crossArea

        print "Layer \t Area (mm2) \t\t Angle (degrees) \t\t Angle (mrad) \t\t R1 (mm) \t\t R5"
        for i in range(len(self.geoArea)):
            print i,"\t",self.crossArea[i],"\t\t",math.degrees(self.geo.getAlpha(self.radiiOfLayers[i][2])),\
                "\t\t",self.geo.getAlpha(self.radiiOfLayers[i][2])*1000,"\t\t",self.radiiOfLayers[i][0],"\t\t",self.radiiOfLayers[i][4]
        #print self.effArea

    def drawOptic(self,focalLength):
        plot = Plot(self.glassLength,self.hDistance,self.mandrelRadius,self.boreDiameter,self.radiiOfLayers,self.d)

    def getAxionSpectrum(self,energy):
        aspectrum = [i**(2.481)*math.exp(-i/1.205) for i in energy]
        maximum = max(aspectrum)
        aspectrum = [2.5*i/maximum for i in aspectrum]
        return aspectrum

    def calculateReflectivityInIDL(self,focalLength):
        self.idl.n = self.n
        self.idl.gamma = self.gamma
        self.idl.dmin = self.slope_d_min
        self.idl.dmax = self.slope_d_max
        self.idl.emin = self.eMin
        self.idl.emax = self.eMax
        self.idl.evalues = self.energyValues
        self.idl.subs = self.subs
        self.idl.mat1R = self.mat2Roughness
        self.idl.mat2R = self.mat2Roughness
        self.idl.subsR = self.subsRoughness
        dE = (self.eMax-self.eMin)/self.energyValues

        self.idl.totalcalc = len(self.n)*len(self.gamma)*len(self.slope_d_min)*len(self.slope_d_max)/2


        recipe_file = open("output/%s/recipe_file.txt" % (self.time), "a")

        print '\n \n'
        print 'Number of calculations: ' + str((self.idl.totalcalc*(len(self.ML_mats[:][0])+len(self.SL_mats)))*len(self.radiiOfLayers))
        print 'Energy resolution: ' + str(self.energyValues) + ' Values over energy range.\n\n'

        self.idl('.Compile LOOPOVERCOATINGS')
        self.idl('.Compile LOOPOVERCOATINGS')
        Aeffective = numpy.zeros((len(self.ML_mats),self.energyValues+1))
        IntRefl = numpy.zeros(len(self.ML_mats))
        fig2 = plt.figure()
        bb = fig2.add_subplot(111)
        for j in range(len(self.ML_mats)):
            print "Materials: " + str(self.ML_mats[j][0]) + "\t" + str(self.ML_mats[j][1])
            self.idl('GETNK,subs,"%s","%s",emin,emax,evalues,NK_ARRAY,LAMBDA,ENERGY,dE' % (str(self.matCombArray[j][0]),str(self.matCombArray[j][1])))
            for i in range(len(self.radiiOfLayers)):
                if i % self.N_per_recipe == 0:

                    alpha = []
                    for k in range(self.N_per_recipe):
                        if (i+k)+1 <= len(self.radiiOfLayers):
                            alpha.append(90-math.degrees(self.geo.getAlpha(self.radiiOfLayers[(i+k)][2])))

                    self.idl.alpha = alpha

                    self.idl('LOOPOVERCOATINGS,n,gamma,dmin,dmax,\
                                alpha,nk_array,lambda,energy,dE,\
                                totalcalc,mat1R,mat2R,subsR,bestml_array,\
                                best_int_array,bestRA,aspectrum,qe,\
                                recipe_refl_array,recipe_ra_array,recipe_int_response_array')



                    for l in range(len(alpha)):
                        recipe_layer = "Layer: " + str(i+l) + ":" + str(self.idl.bestml_array.tolist()) + "\n"
                        print recipe_layer
                        recipe_file.write(recipe_layer)
                        #print self.idl.recipe_refl_array[:,l].tolist()
                        #Integrated_effective = sum([m*dE for m in ])
                        #Integrated_effective = sum([n*dE for n in self.idl.recipe_ra_array[:,l].tolist()])
                        #IntRefl[j] += (Integrated_effective**2)*self.geoArea[i+l]*0.8/10**2
                        ###LayerEffective = [(m**2)*self.crossArea[i+l]*0.8/(10**2) for m in self.idl.recipe_ra_array[:,l].tolist()]
                        LayerEffective = [(m**2)*self.crossArea[i+l]*0.8 for m in self.idl.recipe_ra_array[:,l].tolist()]
                        #print self.idl.recipe_ra_array[10,l].tolist()*self.crossArea[i+l]*0.8/10**2
                        IntRefl[j] += sum([n*dE for n in LayerEffective])
                        Aeffective[j] = [sum(pair) for pair in zip(Aeffective[j],LayerEffective)]
                        a = []
                        #for k in range(len(self.idl.energy)):
                        #    a.append([self.idl.energy[k],self.idl.recipe_ra_array[k,l].tolist()])
                        a = zip(self.idl.energy,self.idl.recipe_ra_array[:,l].tolist())
                        numpy.savetxt('%s/layer%s_%s.txt' % (self.directory+'/'+str(self.ML_mats[j]),i+l,'_'+str(focalLength)),a,delimiter='\t')

                    bb.plot(self.idl.energy,self.idl.best_int_array,
                            label="Layer %s, Int:%s, N:%s, Gamma:%s, dmin:%s, dmax:%s"
                            % (i,'{0:.4g}'.format(self.idl.bestml_array[0]),
                               self.idl.bestml_array[1],self.idl.bestml_array[2],
                               self.idl.bestml_array[3],self.idl.bestml_array[4]))
                    #a = []
                    #for k in range(len(self.idl.energy)):
                     #   a.append([self.idl.energy[k],Aeffective[j][k]])

                    #numpy.savetxt('%s/layer%s_%s.txt' % (self.directory+'/'+str(self.ML_mats[j]),i,'_'+str(focalLength)),a,delimiter='\t')
                    #IntRefl[j] += self.idl.bestml_array[0]

            bb.plot(self.idl.energy,self.idl.aspectrum,label="Axion spectrum")
            bb.plot(self.idl.energy,self.idl.qe,label="Detector quantum efficiency")
            bb.set_xlabel('Energy [keV]')
            bb.set_ylabel('Reflectivity')
            bb.set_title('Material combination: %s, focal length: %s' % (str(self.ML_mats[j]),str(focalLength)))
            bb.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1))
            fig2.savefig('output/%s/%s-%s_%s.png'\
                         % (self.time,str(self.ML_mats[j][0]),str(self.ML_mats[j][1]),str(focalLength)), bbox_inches='tight')
            cla()
            #show()
        return Aeffective,IntRefl

    def plotEffectiveArea(self,MLAeffective):#,SLAeffective):
        for i in range(len(self.ML_mats)):
            plot(self.idl.energy,MLAeffective[i],label='Mat. combo: %s' % (self.ML_mats[i]))
        #for i in range(len(self.SL_mats)):
        #    plot(self.idl.energy,SLAeffective[i],label='SL mat: %s' % (self.SL_mats[i]))

        xlabel('Energy [keV]')
        ylabel('Effective reflectivity * area [mm2]')
        title('Effective reflectivity per area comparison')
        lgd = legend(loc='upper center', bbox_to_anchor=(0.5,-0.1))
        savefig('output/%s/Effective_refl_area.png' % (self.time), bbox_extra_artists=(lgd,), bbox_inches='tight')
        cla()

    def saveEffAreaToFile(self,Eff_area):
        for j in range(len(Eff_area)):
            a = range(len(Eff_area[j]))
            for i in range(len(Eff_area[j])):
                a[i] = [self.idl.energy[i],Eff_area[j][i]]
            numpy.savetxt('output/%s/%s.txt' % (self.time,self.ML_mats[j]),a,delimiter='\t')

    def removeEndOfString(self):
        for i in range(len(self.matCombArray)):
            self.matCombArray[i] = self.matCombArray[i][0][:-10] + '-' + self.matCombArray[i][1][:-10]

        for i in range(len(self.SLMat_array)):
            self.SLMat_array[i] = self.SLMat_array[i][:-10]


    #def plotIntReflComp(self,MLIntRefl,SLIntRefl):
    def plotIntReflComp(self,MLAeffective,SLAeffective,focalLength):

        IntAreaArray = numpy.zeros(len(MLAeffective)+len(SLAeffective))


        for i in range(len(MLAeffective)):
            IntAreaArray[i] = MLAeffective[i]

        for i in range(len(SLAeffective)):
            IntAreaArray[i+len(MLAeffective)] = SLAeffective[i]

        labels = []

        for i in range(len(MLAeffective)):
            labels.append(str(self.ML_mats[0]) + ' ' + str(focalLength[i]) + 'mm')

        for i in range(len(SLAeffective)):
            labels.append(str(self.SL_mats[0]) + ' ' + str(focalLength[i]) + 'mm')

        #NormIntAreaArray = [i/max(IntAreaArray) for i in IntAreaArray]

        #TotalIntRefl = numpy.append(MLIntRefl,SLIntRefl)
        #MaxIntRefl = max(TotalIntRefl)

        #NormTotIntRefl = [i/MaxIntRefl for i in TotalIntRefl]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ind = numpy.arange(len(IntAreaArray))
        width = 0.35

        #rects = ax.bar(ind,NormTotIntRefl,width)
        rects = ax.bar(ind,IntAreaArray,width)
        ax.set_ylim((0,1.1*max(IntAreaArray)))
        ax.set_ylabel('Throughput / sqrt(spot)')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( labels , rotation=-45)
        k= 0
        for i in rects:
                value = IntAreaArray[k]
                height = i.get_height()
                ax.text(i.get_x()+i.get_width()/2., 1.05*height, '{0:.3g}'.format(value),ha='center', va='bottom')
                k += 1

        savefig('output/%s/Relative_eff_throughput.png' % (self.time), bbox_inches='tight')
        #show()


    def calculateSLReflectivityIDL(self,focalLength):
        self.idl.emin = self.eMin
        self.idl.emax = self.eMax
        self.idl.evalues = self.energyValues
        self.idl.D = self.D
        self.idl.subs = self.subs
        self.idl.matR = self.matR
        self.idl.subsR = self.subsRoughness

        dE = (self.eMax-self.eMin)/self.energyValues


        Aeffective = numpy.zeros((len(self.SL_mats),self.energyValues+1))
        IntRefl = numpy.zeros(len(self.SL_mats))

        for i in range(len(self.SL_mats)):
            self.idl('GETSLNK,subs,"%s",emin,emax,evalues,NK_ARRAY,LAMBDA,ENERGY,dE' % (str(self.SLMat_array[i])))
            self.idl.pro('.r SLREFLECTIVITY')
            for j in range(len(self.radiiOfLayers)):
                self.idl('SLREFLECTIVITY,%s,D,matR,subsR,nk_array,lambda,energy,dE,RA,aspectrum,qe,int_array,int_response' % (90-math.degrees(self.geo.getAlpha(self.radiiOfLayers[j][2]))))

                #Integrated_effective = sum([n*dE for n in self.idl.int_array.tolist()])
                #IntRefl[i] += (Integrated_effective*self.geoArea[j])/(10**6)

                #LayerEffective = [k*self.geoArea[j] for k in self.idl.int_array.tolist()]

                #Aeffective[i] = [sum(pair) for pair in zip(Aeffective[i],LayerEffective)]

                LayerEffective = [(k**2)*self.crossArea[j]*0.8/10**2 for k in self.idl.RA ]

                Aeffective[i] = [sum(pair) for pair in zip(Aeffective[i],LayerEffective)]

                #Integrated_effective = sum([n*dE for n in self.idl.RA.tolist()])
                #IntRefl[i] += (Integrated_effective**2)*self.geoArea[j]*0.8/10**2
                IntRefl[i] += sum([n*dE for n in LayerEffective])


        return Aeffective,IntRefl


    def initializeIDL(self):
        self.idl('getnk,"%s","%s","%s",%s,%s,%s,%s,%s,%s,%s'\
                 % (str(self.subs),str(self.mat1),str(self.mat2),self.eMin,
                    self.eMax,self.energyValues,'nk_array','LAMBDA','energy','dE'))


    def makePlotsForLOI(self,focalLength):
        dE = (self.eMax-self.eMin)/(self.energyValues+1)

        energy = numpy.arange(self.eMin,self.eMax,dE)


        Aeffective = numpy.zeros((len(self.ML_mats),self.energyValues+1))
        Aeffective_w_det = numpy.zeros((len(self.ML_mats),self.energyValues+1))
        OpeningRefl = numpy.zeros((len(self.ML_mats),self.energyValues+1))
        OpeningRefl_w_det = numpy.zeros((len(self.ML_mats),self.energyValues+1))



        IntRefl = numpy.zeros(len(self.ML_mats))

        qe = numpy.loadtxt('Eff_mM_No_Strongback.txt', delimiter = '\t')

        det_qe = numpy.interp(energy,qe[:,0],qe[:,1])

        for i in range(len(self.ML_mats)):
            totalcrossarea = 0
            for j in range(len(self.radiiOfLayers))[1:]:

                data = numpy.loadtxt('output/160413-1510/W-B4C/layer%s__%s.txt' % (j,str(focalLength)), delimiter = '\t')

                refl = data[:,1]
                refl_w_detector = data[:,1]*det_qe

                LayerOpeningRefl = [(m**2)*self.crossArea[j]*0.8/(10**2) for m in refl]
                LayerOpeningRefl_w_det = [(m**2)*self.crossArea[j]*0.8/(10**2) for m in refl_w_detector]
                OpeningRefl[i] = [sum(pair) for pair in zip(OpeningRefl[i],LayerOpeningRefl)]
                OpeningRefl_w_det[i] = [sum(pair) for pair in zip(OpeningRefl_w_det[i],LayerOpeningRefl_w_det)]

                totalcrossarea += self.crossArea[j]

                LayerEffective = [(k**2)*self.crossArea[j]*0.8/10**2 for k in refl]
                Aeffective[i] = [sum(pair) for pair in zip(Aeffective[i],LayerEffective)]
                LayerEffective = [(k**2)*self.crossArea[j]*0.8/10**2 for k in refl_w_detector]
                Aeffective_w_det[i] = [sum(pair) for pair in zip(Aeffective[i],LayerEffective)]



        return OpeningRefl[0,:], OpeningRefl_w_det[0,:], Aeffective[0,:], Aeffective_w_det[0,:],totalcrossarea/10**6
if __name__ == '__main__':
    Main()
