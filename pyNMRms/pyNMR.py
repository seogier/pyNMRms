# -*- coding: utf-8 -*-
r"""
Created on Jan 2017
NMR Data processing for NIST MRI Measurement service
Analyzes standard pulse sequences using Tecmag console and file format. Can read in limited data from Bruker and Varian 
Need to install the following packages from command line
    Anaconda3       #this includes python,  scipy, numpy, and lots of other useful packages;  download and run the anaconda installer for the system that you are using
    pyqtgraph            this includes all of the plotting packages used, use the latest version, current version is 
    lmfit                this includes the nonlinear least squares fitting
    PyOpenGL             used for 3d plotting, not critical
    nmrglue              opens Varian and Bruker data
    reportlab            generates pdf reports
    pydicom==1.4.2       for importing and exporitng dicom files, not critical, only package at present that needs old version
 
 CAn run the following lines in anaconda3\scripts where pip is located   
    pip install pyqtgraph           
    pip install lmfit               
    pip install PyOpenGL            
    pip install nmrglue             
    pip install reportlab            
    pip install pydicom==1.4.2        
     
    Current package versions
    Python=3.12.4 | packaged by Anaconda, Inc. | (main, Jun 18 2024, 15:03:56) [MSC v.1929 64 bit (AMD64)];
    PyQt=5.15.10;
    pyqtgraph=0.13.7
    
Uses pyMRIGui created from pyMRIGui.ui by PyQt5.uic.pyuic called from pyuic5.bat
  execute from system shell to regenerate GUIs (note pyuic5.bat needs a correct path to python, which may vary)
    designer\pyuic5 designer\pyNMRGUI.ui -o pyNMRms\pyNMRGUI5.py 
    designer\pyuic5MRLab designer\pyNMRGUI.ui -o pyNMRms\pyNMRGUI5.py 
    
@author: stephen russek

All times in s, magnetic fields in tesla, frequencies in Hz, magnetic gradients in T/m, gradient calibration in T/m/A

filename keys and conventions
    'IR'                       :if found sets dataType to IR
    'CPMG' or 'CPFH'                    :if found sets dataType to CPMG
    'CPMG_tau_array'            :if found sets dataType to CPMG_tau_array
    '_SE_' or 'SpinEcho'       :if found sets dataType to SE
    't1rho'                    :if found sets  dataType='T1rho', case insensitive
    'nutation'                 :if found sets  dataType='nutation', case insensitive
    'PGSE' or 'BPP-LED'        :if found sets dataType to Diffusion
    'PD'                       :if found sets dataType to PD, performs proton density analysis
    '_Gx_' or 'GxCal'            :if found sets gradOrientation='X'    
    '_Gy_' or 'GyCal'            :if found sets gradOrientation='Y'   
    '_Gz_' or 'GzCal'            :if found sets gradOrientation='Z'   
    
Uses a modified version of pytnt to open and extract data from TEcmag tnt files, not fully implemented and does not fully unpack the pulse sequence structure
TNMR parameter tables searched for in module processTNT ((de|lp)[0-9]+:[0-9])|TI_[1-9]|ti_times|cpmgloop|gr0:2|GxTable|GyTable|gztable|rdArray|seloop|HPMW|TX3:2   
    seloop        Loop table for spin echo delay, SE=2*seloop+t180

To make a Windows executable, pip install pyinstaller, run command line below with correct path to pyNMRms.py
 pyinstaller --noconfirm --onedir --windowed --icon "C:/Users/serus/Desktop/pyNMRms/icons/Professor.ico" --add-data "C:/Users/serus/Desktop/pyNMRms/pyNMRms/NISTlogo.jpg;."  "C:/Users/serus/Desktop/pyNMRms/pyNMRms/pyNMR.py"
  or 
 pip install auto-py-to-exe
 run auto-py-to-exe.exe with the switches above selected
"""

VersionID='pyNMRms V9-07-2024'
import sys, os, time, datetime, copy        #useful system libraries
import numpy as np
from scipy import constants
from scipy import integrate, interpolate
#from scipy.integrate import odeint
from scipy.signal import find_peaks_cwt
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from scipy.stats import chisquare
from scipy.stats import f as fdis
import re
import pprint
from pyqt import *         #imports required PyQt modules, tries PyQT4 then PyQt5
import pyqtgraph as pg    #uses pyqtgraph PlotWidget for graphs and ImageView windows for images, http://www.pyqtgraph.org/documentation/introduction.html
import pyqtgraph.opengl as gl
import pyqtgraph.exporters

from numpy import fft
import lmfit
import T1IR, T2CPMG, T2CPMGbiExp, DiffusionPGSE, dampedSin, nutatefit, LorentzComplex, expFit #fitting modules for NMR data based on  lmfit
import image1dfit       #used to fit id images to calibrate gradients
import nmrglue as ng        #nmrGlue is used to read a variety of NMR files including Bruker and Varian
import pydicom    #pydicom is used to import/export DICOM images  pydicom.UID
from pydicom.dataset import Dataset, FileDataset
from TNMRviewer import TNMRviewer
# maintain backward compatibility with QT4
if pyqtVersion==4:
  from pyNMRGUI4 import Ui_pyNMRGui    #main window Gui
if pyqtVersion==5:
  from pyNMRGUI5 import Ui_pyNMRGui    #main window Gui
from uncertainty import UncertaintyWindow
from processTNT import TNTfile
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
#import pyperclip as clipboard  #copy and paste from system clipboard
#from Pillow import Image
#from PIL import Image           #used for resizing 2d images
from scipy.ndimage import zoom   #used for interpolation of images
import imageScaleDialog     #execute: designer\pyuic5 designer\imageScaleDialog.ui -o pyNMRms\imageScaleDialog.py
QApplication.setAttribute(Qt.AA_DisableHighDpiScaling, False) #enable highdpi scaling
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, False) #use highdpi icons





class pyNMR(QMainWindow):
  def __init__(self ,parent = None):
        ''' Class for main pyNMR window (using pyqt gui), designed to import and analyze NMR data, mostly from techmag formated data'''
        super(pyNMR, self).__init__()
        self.ui = Ui_pyNMRGui()
        self.ui.setupUi(self)
        self.setAttribute(Qt.WA_NativeWindow, True) 
        self.bblabelStyle = {'color':'b', 'font-size': '18px'}
        self.bbtitleStyle = {'color':'b', 'font-size': '18px'}
        self.wblabelStyle = {'color':'b', 'font-size': '18px'}
        self.wbtitleStyle = {'color':'b', 'font-size': '18px'}
        self.labelStyle=self.bblabelStyle
        self.titleStyle=self.bbtitleStyle
        self.setWindowTitle(VersionID)
        self.resize(1200,1000)
        self.dataReal=self.ui.gvReal     #pyqtgraph widget for pulse sequence plot
        self.dataReal.showAxis('right')
        self.dataReal.showAxis('top')
        self.dataImg=self.ui.gvImag
        self.dataImg.showAxis('right')
        self.dataImg.showAxis('top')
        self.showMaximized() 
        self.phaseUnwrap=True       #flag to unwrap phase data
        self.ncurveMax=50   #maximum curves to plot, otherwise just plot individual curves so the program does not crash
        self.RIorAP='RI'       #flag to plot Real/Imaginary or Amplitude/Phase
        self.inf1 = pg.InfiniteLine(movable=True, angle=90, label='x={value:0.4e}', 
                labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
        self.inf2 = pg.InfiniteLine(movable=True, angle=0, pen=(0, 0, 200),  hoverPen=(0,200,0), label='y={value:0.4e}', 
                labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
        self.inf3 = pg.InfiniteLine(movable=True, angle=90, label='x={value:0.4e}', 
                labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
        self.inf4 = pg.InfiniteLine(movable=True, angle=0, pen=(0, 0, 200),  hoverPen=(0,200,0), label='y={value:0.4e}', 
                labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
        self.t90=13.5e-6       #pulse duration for 90 degree rotation (s)
        self.ui.let90.setText("{:.4f}".format(self.t90*1E6))
        self.ui.leT180.setText("{:.4f}".format(2E6*self.t90))
        self.dataType='FID'      #dataType can be FID, IR, nutation, CPMG, SE, Diffusion, PD, T1rho, HP RFattn
        self.rawDataType='rawMag'     #Specifies raw data format to be displayed in plots or DICOM export, can be rawMag, rawPhase, reconMag, reconPhase
        self.TIp5s=np.fromstring('0.001s  0.0013869s  0.0019235s  0.0026678s  0.0037s  0.0051316s  0.0071172s  0.0098709s  0.01369s  0.018987s  0.026334s  0.036523s  0.050654s  0.070253s  0.097435s  0.13513s  0.18742s  0.25994s  0.36051s  0.5s',sep='s')
        self.TI1s=np.fromstring('0.001s 0.0014384s 0.0020691s 0.0029764s 0.0042813s 0.0061585s 0.0088587s 0.012743s 0.01833s 0.026367s 0.037927s 0.054556s 0.078476s 0.11288s 0.16238s 0.23357s 0.33598s 0.48329s 0.69519s 1.0s',sep='s')
        self.TI=self.TI1s   #change
        self.TIString=re.compile('IR_[0-9]*s') 
        #self.TIString2='de7:2'
        self.TIString2='ti_times'
        self.T1reported=0.0   #reported T1 time in s
        self.T2reported=0.0   #reported T2 time in s
        self.T1Uncertainty=0.0   #reported T1 uncertainty in s
        self.T2Uncertainty=0.0   #reported T2 uncertainty in s
        self.dT1dT=0.03   #normalized T1 temperature derivative 1/C
        self.dT2dT=0.03   #normalized T2 temperature derivative 1/C
        self.T1FitOK=False    #flags that all need to be true to accept data
        self.T2FitOK=False
        self.T1residualsOK=False
        self.T1recoveryTimeOK=False
        self.T1fwhmTimeOK=False
        self.T2residualsOK=False
        self.T2recoveryTimeOK=False
        self.T2fwhmTimeOK=False
        self.temperatureStable=False
        self.allQCtests='Pass'
        self.acqDelay=5E-6    #standard time delay to open acquisition gate before acquisition
        self.nT1sRequired=5   #number of T1 required to make sure Mz has fully relaxed
        self.recoveryTime=0 #acquisition time + last delay
        self.deltaRequired=0.97 #Inversion efficiency must be above self.deltaRequired
        self.tSEarrayname='seloop'       #delay table for spin echo times
        self.tSEarraynameOld='de2:2'       #delay table for spin echo times        
        #CPMG parameters
        self.cpmgLoopTable=np.fromstring('7  15  23  31  39  47  55  63  71  79  87  95  103  111  119  127  135  143  151  159', sep=' ')
        self.cpmgLoopIt=self.cpmgLoopTable+1
        self.cpmgLoopID='lp0:2'
        self.cpmgLoopID2='cpmgloop'
        self.tauCPMG=0.001      #CPMG delay in s
        self.cpmgT180=2*self.t90        #CPMG 180 pulse
        self.cpmgDelay=5.0E-6           #CPMG RX blank delay 
        self.HPMWTableID='HPMW'     #table to indicate if microwaves are on (0) or microwaves are off
        self.HPMWTableID2='TX3:2'     #table to indicate if microwaves are on (0) or microwaves are off
        #***Pulsed Gradient Spin Echo (PGSE) for diffusion       
        self.bcalfactor=10.0             #converts dac units to bvalues in s/mm^s
        self.gradAmpArray=np.fromstring('0  0.5  1  1.5  2  2.5  3  3.5  4  5  7  10 ', sep =' ')           #array of output DAC value multiplyers, divide by Gmax
        self.PGSEgrad=0      #gradient pulse width in s
        self.ui.leGradPW.setText(str(self.PGSEgrad*1000))
        self.gradPulseRiseTime=0.1E-3       #gradient pulse rise-time in s
        self.ui.leGradRisetime.setText(str(self.gradPulseRiseTime*1000))
        self.GradPulseType='trap'       #gradient pulse type trap, hSin, Sin2
        self.gradIntF=1.0   #factor determining qmax from gradient pulse duration and amplitude, 1 for square, 2/pi for sin ...
        self.PGSEdelta=0         #delay between two PGSE gradient pulses in s
        self.ui.lePGSEdelta.setText(str(self.PGSEdelta*1000))
        self.deltaString=re.compile('\_\d\dmdelta', re.IGNORECASE)    #regular expression to match delta time delay in filename
        self.deltaString2=re.compile('delta=\d+ms', re.IGNORECASE)
        self.DACmaxString=re.compile('dacmax=\d+', re.IGNORECASE) 
        self.tflipString=re.compile('tflip=\d+\.?\d*', re.IGNORECASE)
        self.gradString=re.compile('\_\d\dmgrad', re.IGNORECASE)    #regular expression to match delta time delay in filename
        self.gradString2=re.compile('grad=\d+\.?\d*ms', re.IGNORECASE)
        self.GrAmpString=re.compile('GrAmp=\d+\.?\d*', re.IGNORECASE) #gradient amp scaling factor locator string
        #spectrum analysis
        self.spectraMinimum=0.02        #sets a minimum threshold (signal/signalMax) for which data in T1, T2, ADC spectra will be analyzed

        #Gradient calibrations
        #gradient direction=X, Gcal(mT/m/A)= 48.680, Imax(A)= 7.076 cal = DifSetup\5_1DImage_PGSE_GxCal_16C_20220816_131946.tnt
        #gradient direction=Y, Gcal(mT/m/A)= 47.266, Imax(A)= 7.246 ,cal== DifSetup\6_1DImage_PGSE_GyCal_16C_20220816_141325.tnt
        self.GrAmp=100.0      #gradient amp scaling factor, standard pulse sequence parameter used in Tecmag diffusion sequences
        self.GxCal=48.680E-3             #gradient calibration in T/m/A 
        self.GyCal=47.266E-3         #Nominal probe value 42.5E-3 , 47.69E-3 from 4.2 mm NMR tube cal            #gradient calibration in T/m/A 
        self.GzCal=54.2E-3             #gradient calibration in T/m/A
        self.GCal=50E-3
        self.Axmax=7.076     #Maximum X gradient amp current when TNMR Amplitude is 100
        self.Aymax=7.246     #Maximum Y gradient amp current when TNMR Amplitude is 100
        self.Azmax=7.094     #Maximum Z gradient amp current when TNMR Amplitude is 100
        self.Amax=self.Axmax    #choose default gradient orientation to be X
        self.cmCal=2.5      #waveforms from current monitors are 2.5A/V in differential mode
        self.RFampgain=61   #RFamp gain in dBm
        self.PGSEintergradTimeDelay=60E-6       #time interval between PGSE gradients=2delta+T180+PGSEintergradTimeDelay
        self.setGradOrientation('X')            
        self.pgseLoopID = 'gr0:2'
        self.BPPLEDLoopID = 'GxTable'
        self.eddyGxCorrection=0.0      #correction to measured diffusion coefficient
        self.eddyGyCorrection=0.0      #correction to measured diffusion coefficient
        self.eddyGzCorrection=0.0      #correction to measured diffusion coefficient
        self.ui.leECCX.setText(str(self.eddyGxCorrection*1000)) #currently only using x,y to measure ADC
        self.ui.leECCY.setText(str(self.eddyGyCorrection*1000))
        self.performGradientCalandECC =True     #flag to perform gradient Calibrations and Eddy Current Corrections
        self.DxReported=0.0     #Reported x diffusion coefficient with eddy current correction in mm^2/s
        self.DyReported=0.0    #Reported y diffusion coefficient with eddy current correction in mm^2/s
        self.DzReported=0.0    #Reported z diffusion coefficient with eddy current correction in mm^2/s
        self.DavReported=0.0    #Reported average diffusion coefficient with eddy current correction in mm^2/s, currently x,y average only
        self.Kreported=0.0      #dimensionless value characterizing curvature in signal vs b-value data
        self.doNotFitDBelowSNR=False        #Flag to not fit diffusion data below max signal/self.SNRminfit
        self.FakeDiffusionParameters='3000, 0.001, 0.0005, 0.02, 0.001' #parameters to generate fake bi-exponential diffusion data (bmax, Da, Db, fb, 1/SNR)
        
         
        self.maxInhLineField=0.25E-6    #maximum standard deviation of  inhomogeneous field in tesla
        self.tempString=re.compile('\_-?\d\d?\.?\d?C')    #regular expression to match temperature in filename _dC _ddC, _dd.dC or negative values
        self.ui.leTauCPMG.setText(str(self.tauCPMG*1000))
        self.repeatIndex=0      #3rd dimension index Tecmag parameter 'TD3', usually 0, 1, 2
        self.nRepeats=1         #number of repeats, size of TD3 dimension
        self.ui.sbRepeat.setMinimum(-1)
        self.ui.sbRepeat.setMaximum(2)
        self.nfitpoints=1000    #number of points in the model reference curve
        self.npws=10            #+-number of peak widths to integrate over

        self.intbwWide=400   #wide integration bandwidth in Hz to get most protons seen in typical MRI voxel
        self.intbwNarrow=20  #Narrow integration bandwidth in Hz to get protons in a single peak
        self.intbw=self.intbwWide      #bandwidth to integrate over in Hz
        self.ui.leIntBW.setText('{:6.2f}'.format(self.intbw))
        self.pwidth=np.array([ 200,500])      #smoothing parameter used in peak finding taucp=1E-3#CPMG time delay in s
        self.nCPG=20    #number 0f repetitions in CPMG sequence
        self.Phase=0.0    #phase shift of spectra in degrees relative to input spectra
        self.lofcriteria=0.003    #criteria for accepting or rejecting model for fit if p-value<lofcriteria fit indicates model is not appropriate
        self.maxResidualSTD = 0.003 #will accept data id standard deviation of residuals is below this value
        self.maxSigAveWater=0.0     #used to calculate proton density
        self.maxSigAveWaterList=[]      #list of water signal reference values
        self.PDList=[]      #list of proton density values
 
        #output report as pdf file using reportLab
        self.reportList=[]    #list of stuff to add 
        self.NISTlogo="NISTlogo.jpg"
        im = Image(self.NISTlogo, 7*inch, 1.00*inch)
        self.reportList.append(im)
        self.reportList.append(Spacer(1,0.2*inch))
        self.reportstyles=getSampleStyleSheet()
        self.reportstyles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        self.reportRealImage='images/reportRealImageZ.jpg'     #filenames for temporary storage of plot images that go into the report, these are typically overwritten when writing a new report
        self.reportImgImage='images/reportImgImageZ.jpg'
        self.reportFitImage='images/reportFitImageZ.jpg'
        self.reportRealImageZ=0 #index of image placed in the report
        self.reportImgImageZ=0
        self.reportFitImageZ=0
        
        #Physical Constants
        self.GammaPMHzperT=constants.physical_constants["proton gyromag. ratio over 2 pi"][0]   #in MHz/T 
        self.GammaPHzperT=1E6*constants.physical_constants["proton gyromag. ratio over 2 pi"][0]   #in Hz/T 
        self.GammaWaterProtonRadperT=constants.physical_constants["proton gyromag. ratio"][0]*(1-constants.physical_constants["proton mag. shielding correction"][0])                                                                                
        self.Gamma=self.GammaWaterProtonRadperT
        self.Gammaf=self.Gamma/2/np.pi
        
        self.magneticField=0.0    #magnetic field in tesla
        self.temperature=20.0    #temperature in C
        self.maxInhLinewidth=2*self.Gammaf*self.maxInhLineField   #maximum allowed FWHM of inhomogenous linewidth
        #self.t90=1E12/(float(self.ui.leRFamplitude.text())*self.Gammaf*4.0)
        self.dxgcal=3.000E-3    #width of gradient calibration structure in m
        self.ppmOffset=0.0      #offset to be applied to frequency array when converted to ppm
        self.blRegionStart=0.9       #region of spectra used to calculate baseline, assumes all data beyond self.blRegionStart and below self.blRegionStart should be zero!!
        self.blRegionStop=0.99# we do not look at last data points in TechMag FIDs because they are arbitrarily set to 0
        self.TechMagEndZeros=25     #number of points that are added to end of FID by TechMag digital filtering
        self.SNRminfit=1000     #SNR limit used in selecting datapoint to fit
 
        #pulse sequence initial parameters

        #signals and slots
        #File
        self.ui.actionSaveData.triggered.connect(self.saveData)
        self.ui.actionOpenDataFile.triggered.connect(self.openDataFile)     #Open NMR data file
        self.ui.actionOpenDataFile_2.triggered.connect(self.openDataFile)    #used to open diffusion files as part of the diffusion analysis sequence
        
        self.ui.actionSave_Messages.triggered.connect(self.saveMessages)
        self.ui.actionClear_Plots.triggered.connect(self.clearPlots)
        self.ui.actionImport.triggered.connect(self.importFile)
        #Data
        self.ui.actionChange_Plot_Background.triggered.connect(self.changeBackground)
        self.ui.actionToggle_Plot_Background.triggered.connect(self.toggleBackground)
        self.ui.actionShow_Crosshairs.triggered.connect(self.showCrossHairs)
        self.ui.actionCalculate_Baselines.triggered.connect(self.calculateBaselines)
        self.ui.actionFit_Lorentzians.triggered.connect(self.fitLorentzians)
        self.ui.actionDelete_current_repeat.triggered.connect(self.deleteCurrentRepeat)
        self.ui.actionReset_data.triggered.connect(self.resetData)
        self.ui.actionInvert_data.triggered.connect(self.invertData)
        self.ui.actionRemove_leading_trailing_points.triggered.connect(self.removeLeadingTrailingPoints)
        self.ui.actionSubtract_reference_curve.triggered.connect(self.subtractReferenceCurve)
        self.ui.actionPlot_RI_or_AP.triggered.connect(self.plotRIorAP)
        self.ui.actionPrint_PD_Data.triggered.connect(self.printPDdata)
        self.ui.actionClear_Report.triggered.connect(self.clearReport)
        self.ui.actionSave_Report.triggered.connect(self.printReport)
        self.ui.actionAdd_plots_to_report.triggered.connect(self.addDataPlottoReport)
        #self.ui.actionSave_Report.triggered.connect(self.saveMessagestoPDF)
        self.ui.actionEnter_TI_list.triggered.connect(self.enterTIList)
        self.ui.actionEnter_CPMG_list.triggered.connect(self.enterCPMGList)
        self.ui.actionMake_T1IR_data.triggered.connect(self.makeFakeT1IRData)
        self.ui.actionMake_Diffusion_data.triggered.connect(self.makeFakeDiffusionData)
        self.ui.actionShowDataWindow.triggered.connect(self.showDataWindow)
        self.ui.action3dViewer.triggered.connect(self.view3d)
        self.ui.actionView2d.triggered.connect(self.view2d)
        self.ui.actionUncertainty_Calculation.triggered.connect(self.showUncertainty)
        self.ui.actionShow_TNT_header.triggered.connect(self.printFullTNTHeader)
        #Calibrations
        self.ui.actionNutation.triggered.connect(self.nutation)
        self.ui.actionOpen_Gradient_Calibration_File.triggered.connect(self.openGradientCalFile)
        self.ui.actionECCGradientRingdown.triggered.connect(self.ECCGradientRingdown)
        self.ui.actionNoise.triggered.connect(self.noiseAnalysis)
        #T1T2
        self.ui.actionT1IR.triggered.connect(self.T1IR)
        self.ui.actionT1IRSpectra.triggered.connect(self.T1IRSpectra)
        self.ui.actionT2CPMG.triggered.connect(self.T2CPMG)
        self.ui.actionT2CPMG_tau_array.triggered.connect(self.T2CPMG_tau_array)
        self.ui.actionT2SE.triggered.connect(self.T2SE)
        self.ui.actionT1rho.triggered.connect(self.T1rho)
        self.ui.actionT2CPMGSpectra.triggered.connect(self.T2Spectra)
        self.ui.actionHP_Analysis.triggered.connect(self.HPAnalysis)
        self.ui.actionT1IR_narrow_wide_bandwidth.triggered.connect(self.T1IRbiBW)
        self.ui.actionT2CPMGbiExp.triggered.connect(self.T2CPMGbiExp)
        self.ui.actionT2CPMGBiExp_Narrow_Wide_Bandwidth.triggered.connect(self.T2CPMGbiExpNWBW)
        #Diffusion
        self.ui.actionDiffusionAnalysis.triggered.connect(self.diffusionAnalysis)   #complete diffusion analysis for x,y gradient directions
        self.ui.actionGradientCalibrationXY.triggered.connect(self.gradCalXY)
        self.ui.actionDiffusion_Analysis_BiExp_NWBW.triggered.connect(self.diffusionAnalysisBiExpNWBW)  #complete diffusion analysis for x,y gradient directions using a biexponential model and for both a narrow and wide integration bandwidth
        self.ui.actionKurtosisAnalysis.triggered.connect(self.kurtosisAnalysis)   #complete diffusion analysis for x,y gradient directions        
        self.ui.actionReset_Grad_and_bvalues.triggered.connect(self.resetGradandBvalues)
        self.ui.actionInput_Gradient_current_traces.triggered.connect(self.inputGradientCurrentTraces)
        self.ui.actionUse_measured_b_values.triggered.connect(self.useMeasuredBvalues)  #toggles between measured and calculated b-values
        self.ui.actionUse_calculated_b_values.triggered.connect(self.useCalculatedBvalues)
        self.ui.actionCalculate_Image_Widths.triggered.connect(self.calculateImageWidths)
        self.ui.actionInputPGSEDiffusion.triggered.connect(self.inputPGSEDiffusion)
        self.ui.actionDiffusion_Summary.triggered.connect(self.diffusionSummary)
        self.ui.actionInput_eddycurrent_correction.triggered.connect(self.inputEddyCurrentCorrection)
        self.ui.actionUpdate_Gcal_and_Imax.triggered.connect(self.updateGcalImax) 
        self.ui.actionDiffusion_Spectra.triggered.connect(self.diffusionSpectra)
        
        #Imaging
        self.ui.actionDisplay_2D_Image.triggered.connect(self.displayImage)       
        
        #push buttons
        self.ui.pbResetData.clicked.connect(self.resetData)
        self.ui.pbIntegrate.clicked.connect(self.integrateData)
        self.ui.pbPlotFIDMax.clicked.connect(self.plotSignalMax)
        self.ui.pbFindIntWidth.clicked.connect(self.findIntWidth)
        self.ui.pbFFT.clicked.connect(self.FFTData)
        self.ui.pbSubtractBaseLines.clicked.connect(self.subtractBaselines)
        self.ui.actionPlot_B1_z.triggered.connect(self.plotB1z)
        self.ui.sbRepeat.valueChanged.connect(self.repeatIndexChange)
        self.ui.hsPhase.valueChanged.connect(self.inputPhaseSlider)
        self.ui.lePhase.editingFinished.connect(self.inputPhaseLineEdit)
        self.ui.pbPhase.clicked.connect(self.setPhase)
        self.ui.pgFreqShift.clicked.connect(self.frequencyShiftData)
        self.ui.sbViewCurveN.valueChanged.connect(self.plotData)
        self.ui.cbDataType.currentIndexChanged.connect(self.dataTypeChange)
        self.ui.cbGradOrientation.currentIndexChanged.connect(self.changeGradOrientation)
        self.coilTaper=5.5
        self.ui.leRFCoilTaper.setText(str(self.coilTaper))
        self.coilSensitivity()
        self.uncertaintyWin=UncertaintyWindow(self)   #module to calculate uncertainties in T1, T2 using MonteCarlo Bloch Solver
        
        self.ypen=pg.mkPen('y', width=3)
        self.rpen=pg.mkPen('r', width=3)
        self.gpen=pg.mkPen('g', width=3)
        self.wpen=pg.mkPen('w', width=3)
        self.kpen=pg.mkPen('k', width=3)
        self.bpen=pg.mkPen(color=(100, 100, 255), width=3)
        self.penstep=0  #flag to indicate pen color should advance
        self.dataXLabel='Time(s)'  #default data are FID vs time
        self.plotBackgroundBlack=True   #default to black background in plots, change to white when printing
        self.setPlotBackground(black=self.plotBackgroundBlack)
        self.symb=['o', 's', 'd', 't', '+','star', 'p', 'h', 't1','t2', 't3', 'o', 's', 'd', 't', '+','star', 'p', 'h', 't1','t2', 't3','o', 's', 'd', 't', '+','star', 'p', 'h', 't1','t2', 't3']
        self.view3DColor = QColor(255, 255 ,255 , alpha=10)
        self.view3DBackground = QColor(155, 155 ,255 , alpha=10)
        self.view3DTransparency = 1.75   #set transparency scaling for 3dview, 1 = transparency set by voxel value
        self.view3Dinvert = False    #flag to invert contrast in 3d image
        self.direc=''   #default directory to open data
        self.fileName=''        #fileame of current NMR data file
        self.ui.chbClearData.setChecked(True)
        self.dataClear=False
        self.message('<b>NIST MRI Biomarker Calibration Service: </b>' + 'software ' + VersionID + '; Date of analysis: ' +datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.message('<b>Python=</b>'+ sys.version + ';<b> PyQt=</b>' +PYQT_VERSION_STR+ ';<b> pyqtgraph=</b>' + pyqtgraph.__version__+ ' ('+ pyqtgraph.__file__+')' + ', <b>lmfit</b>=='+ lmfit.__version__ +'\n')    #print Python version being used
        self.message('<b>Using gyromagnetic ratio =</b> {:.6f} 10^8 rad/s/T,{:.6f} MHz/s/T'.format(self.Gamma/1E8,self.Gammaf/1E6))

#*********************Opening and Manipulating Data Files****************************************
  def openDataFile (self, message=''):
    '''Opens Tecmag tnt files, extracts NMR data as a 4D complex-float array, determines pulse sequence parameters as best it can'''
    if message==''or message==False:        #When method called from even passed parameters are boolean false, not understood
        message='Open Data File'
    else:
        self.message(self.formatText(message,bold=True, color='blue') )
#     self.ui.cbDataType.setCurrentIndex(self.ui.cbDataType.findText('FID (s)'))
#     self.dataXLabel='Time(s)'
#     self.dataReal.setLabel('bottom', self.dataXLabel, **self.labelStyle)
#     self.dataImg.setLabel('bottom', self.dataXLabel, **self.labelStyle)
    f = QFileDialog.getOpenFileName(self,message, self.direc, "Tecmag Files (*.tnt)")        #Note self.FileName is a qString not a string
    if type(f)==tuple:    #passes  Qstring with PyQt4 and a tuple of strings with PyQt5
      if f[0]=='':
          return 'cancel'
      self.fileName=f[0]
    else:
      if not f:
          return 'cancel'
      self.fileName=str(f)      #make sure fileName is a string, not a Qstring
    self.message(self.formatText('*'*100,color='blue'))
    self.ui.lblFileName.setText(self.fileName)
    self.direc, self.baseFileName = os.path.split(str(self.fileName))
    dir_and_file_string = os.path.join(os.path.split(self.direc)[1], self.baseFileName)
    self.message('<b>File= </b>' + dir_and_file_string)
    self.tntfile=TNTfile(self.fileName)   #open tecmag tnt file
    self.message('<b>Comment=</b>'+self.tntfile.TNMRComment+'\n')
    filename=self.fileName.split('/')[-1]       #Strip off directroy info, only search short filename for data type info 
    #Determine data type
    if filename.find('IR') >=0:    #use the filename to try to determine the type of data
        self.dataType='IR'
    if filename.lower().find('t1rho') >=0:    #use the filename to try to determine the type of data
        self.dataType='T1rho'
    if filename.lower().find('fid90ecc') >=0:
        self.dataType='ECCringdown'
        try:        # find ringdowntimearray
                self.gradRingdownDelay=self.tntfile.DELAY['gradRingdownDelay']
                self.message('<b>gradRIngdownDelay table (s): </b>' + str(self.gradRingdownDelay))
        except:
            self.message('<b>Cannot find gradRingdownDelay table (s): </b>')
    if filename.lower().find('nutation') >=0:
        self.dataType='nutation'
    if filename.find('CPMG') >=0 or filename.find('CPFH') >=0:
        self.dataType='CPMG'
        self.tauCPMG=self.tntfile.tau*1E-3       #The table offset is usually set to 100ns
        self.ui.leTauCPMG.setText(str(self.tauCPMG*1000))
        self.message('<b>CPMG tau(ms) set to: </b>' + str(self.tauCPMG*1000))
    if filename.find('CPMG_tau_array') >=0 :
        self.dataType='CPMG_tau_array'
        self.tauCPMG=self.tntfile.tau*1E-3       #The table offset is usually set to 100ns
        self.ui.leTauCPMG.setText(str(self.tauCPMG*1000))
        self.message('<b>CPMG tau array offset(ms) set to: </b>' + str(self.tauCPMG*1000))
    if filename.find('_SE_') >=0 or filename.find('SpinEcho')>=0  or filename.find('T2SE_NMR')>=0:        #setup for spinecho analysis
        self.dataType='SE'
        if filename.find('T2SE_NMR')>=0:
            self.tSEarrayname='teDelay' #the T1SE_NMR sequence usies a different SE array table anme to be consistent with other MRI sequences
        self.tauSEstep=0.0      #time step for echo time array
        ind=filename.find('_tau=')
        if ind >=0:
            self.tauSEstep= float(re.search('tau=(.+?)ms', filename).group(1))/1000 #spin echo times in sec
        self.tauSEstart=0.0      #initil echo time for echo time array
        ind=filename.find('_tauStart=')
        if ind >=0:
            self.tauSEstart= float(re.search('tauStart=(.+?)ms', filename).group(1))/1000 #spin echo times in sec
    if filename.find('PGSE') >=0:        #Pulsed Gradient Spin Echo
        self.dataType='Diffusion'
        if filename.find('_Gx_') >=0 or filename.find('GxCal') >=0:
            self.setGradOrientation('X')
        if filename.find('_Gy_') >=0 or filename.find('GyCal') >=0:
            self.setGradOrientation('Y')
        if filename.find('_Gz_') >=0 or filename.find('GzCal') >=0:
            self.setGradOrientation('Z')
        if filename.find('_hSin') >=0: #hSin flag indicates gradient pulses are half sin waves, default is trapezoidal
            self.GradPulseType='hSin'
            self.gradPulseRiseTime=0.0
             
        else:
            self.GradPulseType='trap'
            self.gradIntF=1
        try:        # try to find delta and grad in filename , obsolete, now pulling from tecmag file
#            delta=self.deltaString.findall(filename)[0].replace('mdelta=','').replace('_','')
#            self.PGSEdelta=float(delta)*1E-3     #gradient delay pulse width in s
            self.PGSEdelta=self.tntfile.Delta/1000
            self.ui.lePGSEdelta.setText(str(self.PGSEdelta*1000))
#            grad=self.gradString.findall(filename)[0].replace('_','').replace('mGrad','')
#            self.PGSEgrad=float(grad)*1E-3     #gradient pulse width in s
            self.PGSEgrad=self.tntfile.Grad/1000
            self.ui.leGradPW.setText(str(self.PGSEgrad*1000))
        except:
            try:
                delta=self.deltaString2.findall(filename)[0].replace('delta=','').replace('ms','')
                self.PGSEdelta=float(delta)*1E-3     #gradient delay pulse width in s
                self.ui.lePGSEdelta.setText(str(self.PGSEdelta*1000))
                grad=self.gradString2.findall(filename)[0].replace('grad=','').replace('ms','')
                self.PGSEgrad=float(grad)*1E-3     #gradient pulse width in s
                self.ui.leGradPW.setText(str(self.PGSEgrad*1000))  
            except:
                self.PGSEdelta=0
                self.PGSEgrad=0
                self.message('Did not find PGSE parameters')
    if filename.find('BPP-LED') >=0:     #bipolar pulse-low eddy currents
        self.dataType='Diffusion'
    if filename.find('_PD_') >=0:     #bipolar pulse-low eddy currents
        self.dataType='PD'
    index = self.ui.cbPulseSequence.findText(self.dataType, Qt.MatchFixedString)     #set data type in gui combo box
    self.message('Data type= ' + self.dataType)
    if index >= 0:
        self.ui.cbPulseSequence.setCurrentIndex(index)
    #find useful parameters embeded in the filename
    if self.tempString.search(filename):    #use the filename to try to determine temperature
        t=self.tempString.findall(filename)[0].replace('p','.').replace('C', '').strip('_')
        self.temperature=float(t)
        self.ui.leTemperature.setText(t)
    if self.TIString.search(filename):    #use the filename to try to determine TI list
        self.TIDelayTable=str(self.TIString.findall(filename)[0].strip('s').replace('IR', 'TI')) 
    self.ui.lblF0.setText("{:.6f}".format(self.tntfile.ob_freq[0]))     #"{:.6e}".format(
    self.magneticField=1E6*self.tntfile.ob_freq[0]/self.Gammaf
    self.ui.leField.setText("{:.6f}".format(self.magneticField))
    self.message('<b>Field(T)=</b>' + self.ui.leField.text() +', <b>Obs. Frequency(MHz)=</b>' + self.ui.lblF0.text()+', <b>Temperature(C)=</b>' + self.ui.leTemperature.text())    
#    self.message('<b>t90</b>(&mu;s)=' + self.ui.let90.text() +', <b>t180</b>(&mu;s)=' + self.ui.leT180.text(), report=True)
    self.message('<b>Data Acquisition time (UTC):</b> start= ' + self.tntfile.start_time.isoformat()+ 
                 ', finish= ' + self.tntfile.finish_time.isoformat()) #+', acq time TR (s)= '+ str(self.tntfile.spec_acq_time()))
    #self.datafile = ndm.nmrData(fileName, "ntnmr")
    self.data=self.tntfile.DATA   #input 4 dimensional data array, usually 4th dimension not use
    if filename.find('T1IR_NMR')>=0 or filename.find('T2SE_NMR')>=0:     #MRI sequence with parmagter on 4th dimension, switch it to 2nd dimension
        self.data=np.swapaxes(self.data,1,3)
        self.message('MRI sequence, swapping axes 1 and 3', color='red')
    self.tfid=self.tntfile.fid_times()      #time points for FID waveforms
    self.freq=self.tntfile.freq_Hz()        #frequency points for spectra
    QApplication.processEvents()
    self.message('<b>Data shape:</b>' + ' Number of data points=' + str(self.data.shape[0])+ ', parameters=' + str(self.data.shape[1])+ ', repeats=' + str(self.data.shape[2]))
    self.nPoints=self.data.shape[0]
    self.ui.lenPoints.setText(str(self.nPoints))
    self.nSpectra=self.data.shape[1]
    self.ui.lenSpectra.setText(str(self.nSpectra))
    self.ui.sbViewCurveN.setMaximum(self.nSpectra-1)
    self.nRepeats=self.data.shape[2]
    self.ui.lenTD3.setText(str(self.nRepeats))    
    self.ui.sbRepeat.setMaximum(self.nRepeats-1)
    self.ui.leiStart.setText('0')
    self.ui.lejStop.setText(str(self.nPoints))
    #self.data=self.datafile.allFid[0]
    self.message('<b>T90(&mu;s)</b>= ' + "{:.2f}".format(self.tntfile.T90) + '; <b>T180(&mu;s)</b>= ' + "{:.2f}".format(self.tntfile.T180)+ '; <b>tau(ms)</b>= ' + "{:.2f}".format(self.tntfile.tau))
    self.ui.let90.setText("{:.2f}".format(self.tntfile.T90))
    self.ui.leT180.setText("{:.2f}".format(self.tntfile.T180))
    self.ui.leTauCPMG.setText("{:.2f}".format(self.tntfile.tau))
    self.ui.leNutationIncrement.setText("{:.2f}".format(self.tntfile.NutIncrement))
    self.GrAmp=self.tntfile.GrAmp
    self.message('<b>GrAmp</b>= ' + "{:.2f}".format(self.GrAmp) + ' : Igr=Imax*GrAmp*GrDAC/1E4')
    self.Phase=0.0    #Overal phase shift of entire dataset
    self.ui.hsPhase.setValue(0)   #reset data phase to 0
    self.ui.lePhase.setText('0.0')
    self.PhaseArray=np.zeros((self.nSpectra,self.nRepeats))   #phase shift for specific spectra
    self.FreqArray=np.zeros((self.nSpectra,self.nRepeats))   #frequency of maximum peak for each spectra

    #look in the tnt file delay table dictionary to find useful arrays
    try:        # find TI times BossWay
        self.TI=self.tntfile.DELAY[self.TIDelayTable]+self.acqDelay
        self.message('<b>TI table + Acq  delay (s): </b>' + str(self.TI))
    except:
        #self.message('***Cannot find TI table in Filename***')
        try:        #KarlWay
            self.TI=self.tntfile.DELAY[self.TIString2]+self.acqDelay
            self.message('<b>TI table + Acq  delay (s): </b>' + str(self.TI))
        except:
            if self.dataType=='IR':
                self.message('***Cannot find TI table using standard TI strings or in filename***')
    try:
        self.cpmgLoopTable=self.tntfile.DELAY[self.cpmgLoopID]
        self.cpmgLoopIt=self.cpmgLoopTable+1        #add 1 since the loop table 0 means one interation
        self.message('<b>CPMG loop table:</b> ' + str(self.cpmgLoopIt))
    except:
        pass
    try:
        self.cpmg_tau_array=self.tntfile.DELAY['tau_array']
    except:
        pass
    try:        #try to find the SE loop table
        self.tauSE=self.tntfile.DELAY[self.tSEarrayname]     #tnt arrays can be interpreted as floats
        if self.tauSEstep!=0:
            self.tauSE=(np.arange(self.nSpectra))*self.tauSEstep+self.tauSEstart      #for case wheter SE is stepped and not taken from a table
        self.message('<b>Spin Echo delays:</b> ' + str(self.tauSE))
    except:
        pass
    try:        #try to find the CPMG loop table
        self.cpmgLoopTable=self.tntfile.DELAY[self.cpmgLoopID2].astype(int)     #tnt arrays can be interpreted as floats
        self.cpmgLoopIt=self.cpmgLoopTable+1
        self.message('<b>CPMG loop table:</b> ' + str(self.cpmgLoopIt))
    except:
        pass
    try:        #try to find the Hyperpolarization MW table
        self.HPMWTable=self.tntfile.DELAY[self.HPMWTableID].astype(int)     #tnt arrays can be interpreted as floats
        hpextend= int(self.nSpectra/self.HPMWTable.size)
        self.HPMWmask=np.tile(self.HPMWTable, hpextend)
        self.message('<b>Hyperpolarization table found:</b> ' + str(self.HPMWTable) + ': Extended *' + str(hpextend))
        self.dataType='HP'
        index = self.ui.cbPulseSequence.findText(self.dataType, Qt.MatchFixedString)     #set data type in gui combo box
        if index >= 0:
          self.ui.cbPulseSequence.setCurrentIndex(index)
        
    except:
        pass
    try:        #try to find the Hyperpolarization MW table
        self.HPMWTable=self.tntfile.DELAY[self.HPMWTableID2].astype(int)     #tnt arrays can be interpreted as floats
        hpextend= int(self.nSpectra/self.HPMWTable.size)
        self.HPMWmask=np.tile(self.HPMWTable, hpextend)
        self.message('<b>Hyperpolarization table found:</b> ' + str(self.HPMWTable) + ': Extended *' + str(hpextend))
        self.dataType='HP'
        index = self.ui.cbPulseSequence.findText(self.dataType, Qt.MatchFixedString)     #set data type in gui combo box
        if index >= 0:
          self.ui.cbPulseSequence.setCurrentIndex(index)
        
    except:
        pass
#    try:        #try to find the PGSE loop table
    if filename.find('PGSE') >=0:        #Pulsed Gradient Spin Echo
        self.gradAmpArray=self.tntfile.DELAY[self.pgseLoopID].astype(float)     #tnt arrays can be interpreted as floats
        self.message('<b>PGSE grad amplitude table:</b> ' + str(self.gradAmpArray) +', <b>PGSE gradient pulse type: </b> ' + self.GradPulseType)
        self.gradientCurrent=self.GrAmp/100*self.gradAmpArray/100*self.Amax
        self.gradientAmplitude=self.gradientCurrent*self.GCal   #in T/m
        self.message( '<b>gradient strength (mT/m)=</b>' + np.array2string(self.gradientAmplitude*1000, precision=2))
        if self.PGSEdelta>0 and self.PGSEgrad>0:        #use standard ST formula to calulate b-values
            self.gradPulseRiseTime=self.tntfile.gradPulseRiseTime
            if self.gradPulseRiseTime<0:       #0 indicates that the risetime could not be read in from theTMR file so use the value listed in the user settable box
                if self.GradPulseType=='hSin':
                    self.gradPulseRiseTime=0
                else:
                    self.gradPulseRiseTime=float(self.ui.leGradRisetime.text())/1000
            self.ui.leGradRisetime.setText(str(self.gradPulseRiseTime*1000))
            gPulseWidth=self.PGSEgrad+2*self.gradPulseRiseTime
            gPulseDelay=gPulseWidth+2*self.PGSEdelta+2*self.t90+2*self.PGSEintergradTimeDelay
            self.bValueArray=self.STbvalue(g=self.gradientAmplitude, delta=gPulseWidth, Delta=gPulseDelay, risetime=self.gradPulseRiseTime, pulsetype=self.GradPulseType)
            qmax=self.gradIntF*gPulseWidth*np.amax(self.gradientAmplitude)*self.Gamma/1000     # appoximate maximum q value 
            self.message('<b>qmax (1/mm)=</b> {:.3f}'.format(qmax))
            #self.bValueArray=1E-6*(self.Gamma*self.gradientAmplitude*self.PGSEgrad)**2*(2.0*self.PGSEgrad/3+2*self.PGSEdelta+2*self.t90+2*self.PGSEintergradTimeDelay)
        else:
            self.bValueArray=self.gradAmpArray**2
                #b-value in s/mm^2   so ADC in mm^2/s!!!!
        self.message( '<b>b-values (s/mm^2)=</b> ' + np.array2string(self.bValueArray, precision=2, separator=','))         #"{:.4f}".format(

    if filename.find('BPP-LED') >=0:     #bipolar pulse-low eddy currents
        self.gradAmpArray=self.tntfile.DELAY[self.BPPLEDLoopID].astype(float)     #tnt arrays can be interpreted as floats
        self.message('<b>BPP-LED loop table:</b> ' + str(self.gradAmpArray))
        gradientAmplitude=self.gradAmpArray*self.Amax*self.GCal/100.0   #in T/m
        #self.bValueArray=1E-6*(self.Gamma*gradientAmplitude*self.PGSEgrad)**2*(self.PGSEdelta)
        self.bValueArray=self.gradAmpArray**2        #needs appropriate calibration
#    except:
#        pass
    self.dwellTime=self.tntfile.dwell[0]
    self.lastDelay=self.tntfile.last_delay
    self.ui.leLastDelay.setText("{:.3f}".format(self.lastDelay))
    self.recoveryTime=self.lastDelay+self.dwellTime*self.nPoints
    self.message('<b>Dwell time(s)=</b>' + str(self.dwellTime) + ', <b>Recovery time(s)=</b>' + str(self.recoveryTime))
#    calcualte maximum FID signal used in proton density claculation
    maxSig='Max Signal:'
    self.maxSigAve=0.0
    for j in range(self.nRepeats):
        for i in range(self.nSpectra):
            ms = float(np.amax(np.abs(self.data[:,i,j,0])))
            self.maxSigAve+= ms
            maxSig+= str(j) + ', ' + str(i) +'=' + "{:.3e}".format(ms)
    #self.message(maxSig)
    self.maxSigAve=self.maxSigAve/(self.nRepeats*self.nSpectra)

    if self.dataType=='PD':     #for proton density data record the maximum FID signal and normalize to water signal to get PD
        if filename.find('Water') >=0:
            self.maxSigAveWater=self.maxSigAve
            self.maxSigAveWaterList.append(self.maxSigAve)
        else:
            if self.maxSigAveWater >0:
                self.PD=100*self.maxSigAve / self.maxSigAveWater
                self.PDList.append(self.PD)
                self.message('<b>PD(%)=</b>' + "{:.3f}".format(self.PD))
    QApplication.processEvents()
    self.plotData()
  
  def importFile (self, message=''):
    '''Opens jdx, Bruker files, extracts NMR data as a 4D complex-float array, determines pulse sequence parameters as best it can'''
    if message==''or message==False:        #When method called from even passed parameters are boolean false, not understood
        message='Open Data File'
    self.ui.cbDataType.setCurrentIndex(self.ui.cbDataType.findText('FID (s)'))
    self.dataXLabel='Time(s)'
    self.dataReal.setLabel('bottom', self.dataXLabel, **self.labelStyle)
    self.dataImg.setLabel('bottom', self.dataXLabel, **self.labelStyle)
    f = QFileDialog.getOpenFileName(self,message, self.direc, "JCAMP-DX (*.jdx)")        #Note self.FileName is a qString not a string
    if type(f)==tuple:    #passes  Qstring with PyQt4 and a tuple of strings with PyQt5
      if f[0]=='':
          return 'cancel'
      self.fileName=f[0]
    else:
      if not f:
          return 'cancel'
      self.fileName=str(f)      #make sure fileName is a string, not a Qstring
    self.message(self.formatText('*'*100,color='blue'))
    self.ui.lblFileName.setText(self.fileName)
    self.direc, self.baseFileName = os.path.split(str(self.fileName))
    dir_and_file_string = os.path.join(os.path.split(self.direc)[1], self.baseFileName)
    self.message('<b>File= </b>' + dir_and_file_string)
    self.jdxFile= ng.jcampdx.read(self.fileName) 
    print (self.jdxFile[0])
    self.data=self.jdxFile[1]
    print (self.data)
    self.plotData()
       
  def showUncertainty(self):
    self.uncertaintyWin.show()    
    self.uncertaintyWin.ui.leT1Reported.setText(str(self.T1reported))
    self.uncertaintyWin.ui.leT2Reported.setText(str(self.T2reported))
    self.uncertaintyWin.ui.ledT1dT.setText(str(self.dT1dT))
    self.uncertaintyWin.ui.ledT2dT.setText(str(self.dT2dT))

  def clearReport(self):
    self.reportList=[]
    im = Image(self.NISTlogo, 7*inch, 1.00*inch)
    self.reportList.append(im)
    self.reportList.append(Spacer(1,0.2*inch))
    self.reportRealImageZ =0
    self.reportImgImageZ =0
    self.reportFitImageZ =0
    self.ui.txtMessages.clear()
    self.message('<b>NIST MRI Biomarker Calibration Service: </b>' + 'software ' + VersionID + '; Date of analysis: ' +datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
              
  def printReport(self):
#    fn=os.path.splitext(str(self.fileName))[0]   #set default filename to data file name with pdf extension
    fn="MRIMeasurementService_" + self.ui.leSample.text()+datetime.datetime.now().strftime("_%Y-%m-%d")
    f= QFileDialog.getSaveFileName(self,"Write pdf report file", fn, "Report Files (*.pdf)")
    if not f:
        return None
    if type(f)==tuple:    #passes  string with PyQt4 and a tuple with PyQt5
      self.reportFileName=f[0]
      if self.reportFileName=='':   #file open canceled
          return
    else:
      self.reportFileName=f
    fileName=str(self.reportFileName)
    self.Report = SimpleDocTemplate(fileName,pagesize=letter,rightMargin=42,leftMargin=42,topMargin=12,bottomMargin=18)
    #text=str(self.ui.txtMessages.toPlainText())
    reportListout=copy.copy(self.reportList)
    #try:
    self.Report.build(reportListout)
    #except:
    #  QMessageBox.critical(self, "Message", "Can not open report")
       
  def coilSensitivity(self):
        self.coilTaper=float(self.ui.leRFCoilTaper.text())
        self.coilz=np.arange(-25.,25.,0.1)
        b1v=np.vectorize(nutatefit.b1)
        self.b1z=b1v(self.coilz,1,self.coilTaper)     #define RF coil sensitivity function

  def dataTypeChange(self):
      self.dataReal.clear()
      self.dataImg.clear()
      try:
          self.plotData()
      except:
          pass
  
  def repeatIndexChange(self):
      self.repeatIndex=self.ui.sbRepeat.value()
      self.plotData()
           
  def addCrossHairs(self):
      self.inf1.setValue(0.0)
      self.inf2.setValue(0.0)
      self.dataReal.addItem(self.inf1)
      self.dataReal.addItem(self.inf2)
      self.dataImg.addItem(self.inf3)
      self.dataImg.addItem(self.inf4)
      
  def hideCrossHairs(self):
      self.inf1.hide()
      self.inf2.hide()
      self.inf3.hide()
      self.inf4.hide()

  def showCrossHairs(self):
      self.inf1.show()
      self.inf2.show()
      self.inf3.show()
      self.inf4.show()  
                                        
  def plotRIorAP(self):
    '''set flag to plot Real/Imaginary or Amplitude/Phase'''
    if self.RIorAP=='RI':
        self.RIorAP='AP'
    else:
        self.RIorAP='RI'
    self.dataImg.clear()
    self.dataReal.clear()
    try:
        self.plotData() 
    except:
        pass 
      
  def printPDdata(self):
      self.PDarray=np.array(self.PDList)
      self.waterArray=np.array(self.maxSigAveWaterList)
      self.waterArray/=np.mean(self.waterArray)
      self.message("<b>Proton Density summary:</b>")
      self.message('PD(%)=' + np.array2string(self.PDarray,precision=2))
      self.message('Water signal=' + np.array2string(self.waterArray,precision=4))
      self.pw=self.plotWindow(self)
      self.pw.win.setWindowTitle('Proton Density ' + self.fileName)
      self.pw.dplot.setLabel('bottom', "Sample Number")
      self.pw.dplot.setLabel('left', "Proton Density (%)")
      self.pw.show()
      self.pw.dplot.plot(self.PDarray, pen=None, symbol=self.symb[0],symbolSize=14,)
      
                    
  def plotData(self, showCrossHairs=True):
    '''plots NMR data vs index, time, freq,'''  
    if showCrossHairs:
        self.showCrossHairs()
    else:
        self.hideCrossHairs()
    self.iStart=int(self.ui.leiStart.text())
    self.jStop=int(self.ui.lejStop.text())
    self.dataClear=self.ui.chbClearData.isChecked()
    if self.dataClear:
        self.dataImg.clear()
        self.dataReal.clear()
        self.penstep=0
    if self.RIorAP=='RI':
        self.dataReal.setTitle('Real repeat=' +str(self.repeatIndex), **self.titleStyle) 
        self.dataImg.setTitle('Img repeat=' +str(self.repeatIndex), **self.titleStyle)
    if self.RIorAP=='AP':
        self.dataReal.setTitle('Amplitude repeat=' +str(self.repeatIndex), **self.titleStyle) 
        self.dataImg.setTitle('Phase repeat=' +str(self.repeatIndex), **self.titleStyle) 
        magData=np.absolute(self.data[:,:,:,0])
        phaseData=np.angle(self.data[:,:,:,0])  
    npoints=len(self.data[:,0,0,0])
    npoints90=int(0.9*npoints)
    ncurve=self.ui.sbViewCurveN.value()
    nrepeat=self.ui.sbRepeat.value()
    self.SNR=np.zeros((self.nSpectra,self.nRepeats))
    self.baseline=np.zeros((self.nSpectra,self.nRepeats))
    signalmax=np.max(self.data.real)
    if self.ui.cbDataType.currentText()=='FID (s)':
        x=self.tfid
        self.dataXLabel="Time(s)"
    elif self.ui.cbDataType.currentText()=='Spectra (Hz)':
        x=self.freq
        self.dataXLabel="Frequency(Hz)"
    elif self.ui.cbDataType.currentText()=='Spectra (ppm)':
        x=self.freq/self.tntfile.ob_freq[0]-self.ppmOffset     #convert to ppm with 
        self.dataXLabel="Frequency(ppm)"
    else:
        x=np.arange(npoints)
        self.dataXLabel="Index"
    for j in range(self.nRepeats):
        if j==nrepeat or nrepeat==-1:
          for i in range(self.nSpectra):
            if i==ncurve or (ncurve==-1 and i<self.ncurveMax):
              p=pg.mkPen(pg.intColor(i+self.penstep), width=2)
              if self.RIorAP=='RI':
                self.dataReal.plot(x, self.data[:,i,j,0].real, pen=p, width=2, name='real' + str(j))   #plot real
                self.dataImg.plot(x, self.data[:,i,j,0].imag, pen=p, name='img'+ str(j))   #plot imag
              if self.RIorAP=='AP':
                self.dataReal.plot(x, magData[:,i,j], pen=p, width=2, name='Mag' + str(j))   #plot magnitude
                if ncurve!=-1:
                    if self.phaseUnwrap:
                        self.dataImg.plot(x,np.unwrap(phaseData[:,i,j]), pen=p, name='Phase'+ str(j))   #plot phase, but only if a single curve since it is super slow to plot all  
                    else:
                        self.dataImg.plot(x,phaseData[:,i,j], pen=p, name='Phase'+ str(j))   #plot phase, but only if a single curve since it is super slow to plot all  
              self.SNR[i,j]= signalmax/np.std(self.data[int(0.9*npoints):,i,j,0].real)
              self.baseline[i,j]= np.mean(self.data[int(0.9*npoints):,i,j,0].real)
    self.dataReal.setLabel('bottom',self.dataXLabel , **self.labelStyle)
    self.dataReal.setLabel('left', 'Signal', **self.labelStyle)
    self.dataImg.setLabel('bottom', self.dataXLabel, **self.labelStyle)
    self.dataImg.setLabel('left', 'Signal', **self.labelStyle)
    self.addCrossHairs()
    self.penstep+=10 
       
  def plotFFTData(self): 
    self.dataImg.clear()
    self.dataReal.clear()
    self.dataReal.setTitle('Real TD=' +str(self.repeatIndex), **self.titleStyle) 
    self.dataImg.setTitle('Img TD=' +str(self.repeatIndex), **self.titleStyle)
    if self.ui.cbFitAllTD3.isChecked():
      for j in range(self.nRepeats):
        for i in range(self.nSpectra):
          p=pg.mkPen(pg.intColor(i), width=2)
          self.dataReal.plot(self.ftData[:,i,j,0].real, pen=p, width=2, name='real' + str(j))   #plot real
          self.dataImg.plot(self.ftData[:,i,j,0].imag, pen=p, name='img'+ str(j))   #plot imag 
    else:   
      for i in range(self.nSpectra):
        p=pg.mkPen(pg.intColor(i), width=2)
        self.dataReal.plot(self.ftData[:,i,self.repeatIndex,0].real, pen=p, width=2, name='real' + str(self.repeatIndex))   #plot real
        self.dataImg.plot(self.ftData[:,i,self.repeatIndex,0].imag, pen=p, name='img'+ str(self.repeatIndex))   #plot imag 

  def clearPlots(self):
    self.dataImg.clear()
    self.dataReal.clear()
    self.ui.hsPhase.setValue(0)
    self.ui.lePhase.setText('0.0')
    self.Phase=0.0
    self.PhaseArray=0
    self.ui.txtMessages.clear()
    self.message(VersionID + '; ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))    
    
  def inputPhaseLineEdit(self):
    phase=float(self.ui.lePhase.text())
    self.ui.hsPhase.setValue(int(phase*10))
    self.phaseData(phase)
    
  def inputPhaseSlider(self):  
    phase=float(self.ui.hsPhase.value())/10
    self.ui.lePhase.setText(str(phase))
    self.phaseData(phase)
    
  def phaseData(self, phasedeg):
      '''mulitpies data by a phase factor  np.exp(1j*dphase), 
      picks data to be phase shifted by spectra and repeat selection boxes and Phaseall radio button, 
      if these values are -1, all spectra and repeats are phase shifted'''
      
      if self.ui.rbPhaseAll.isChecked:
          dphase=(phasedeg-self.Phase)*np.pi/180.0
          self.data[:,:,:,0]*=np.exp(1j*dphase)
          self.Phase=phasedeg
      else:
          ncurve=self.ui.sbViewCurveN.value()
          nrepeat=self.ui.sbRepeat.value()
          for j in range(self.nRepeats):
              if nrepeat==-1 or j==nrepeat:
                  for i in range(self.nSpectra):
                      if ncurve==-1 or i==ncurve:
                          dphase=(phasedeg-self.PhaseArray[i,j])*np.pi/180.0
                          self.data[:,i,j,0]= self.data[:,i,j,0]*np.exp(1j*dphase)
                          self.PhaseArray[i,j]=phasedeg
      
      self.plotData()  #plot data
      
  def frequencyShiftData(self):
      '''mulitpies data by a frequency shift  np.exp(1j*2pift)'''
      f=float(self.ui.leFreqShift.text()) 
      for j in range(self.nRepeats):
        for i in range(self.nSpectra):
            self.data[:,i,j,0]= self.data[:,i,j,0]*np.exp(1j*2*np.pi*f*self.tfid)
      self.plotData()  #plot data
               
  def integrateData(self):
      '''integrates real part of spectra over desired range to get a single value
      Plots data against proper parameter determined from pulse seqeunce type'''
      self.dataType=self.ui.cbPulseSequence.currentText()
      self.iStart=int(self.ui.leiStart.text())
      self.jStop=int(self.ui.lejStop.text())
      self.repeatIndex=self.ui.sbRepeat.value()
      nspectra=int(self.ui.lenSpectra.text())    #     self.data.shape[1]
      nrepeats=self.data.shape[2]
      self.message('<b>Integrating real part</b> of ' + self.dataType + ' data from ' + '{}'.format(self.iStart) + ' to ' + '{}'.format(self.jStop) + ',  from ' + '{:6.2f}Hz'.format(self.freq[self.jStop]) + ' to ' + '{:6.2f}Hz'.format(self.freq[self.iStart]) + ', <b>Phase adj</b>(degree)=' + self.ui.lePhase.text() + '\n')
      self.dataInt=np.zeros((nspectra,nrepeats))  #nTMR data is 4d, here we are using only 3d data, dim 0 is FID/spectrum, 1 is priamery parameter TI,tcp, 2 is measurement repeats
      #if self.ui.cbFitAllTD3.isChecked():
      if self.dataType=='FID' or self.dataType=='HP':
          x=np.arange(nspectra) 
      if self.dataType=='IR' or  self.dataType=='T1rho':
          if len(self.TI) != nspectra:
              self.message('Data does not match T1 array')
              return
          x=self.TI[:nspectra]
      if self.dataType=='nutation':
          x=np.arange(nspectra) 
      if self.dataType=='SE':
          self.T180=float(self.ui.leT180.text())*1E-6       #need to find teh T180 RF pulse time to add to get the echo time
          self.tSE=2*self.tauSE+self.T180
          x=self.tSE
      if self.dataType=='CPMG':     #for CPMG data convert loop table into aquaition times
          self.tauCPMG=float(self.ui.leTauCPMG.text())*1E-3     #read in tauCPMG in ms convert to s
          self.cpmgT180=float(self.ui.leT180.text())*1E-6
          self.taCPMG=(2.0*(self.tauCPMG+self.cpmgDelay) +self.cpmgT180)*self.cpmgLoopIt
          x=self.taCPMG[:nspectra]
      if self.dataType=='CPMG_tau_array':     #for CPMG data convert loop table into aquaition times
          self.cpmgT180=float(self.ui.leT180.text())*1E-6
          print ('tauarray', self.cpmg_tau_array, 'Offset',self.tauCPMG )
          self.tRefocusCPMG=(2.0*(self.cpmg_tau_array+self.tauCPMG+self.cpmgDelay) +self.cpmgT180)     # Create array of refocus times  *self.cpmgLoopIt
          self.taCPMG=np.outer(self.cpmgLoopIt,self.tRefocusCPMG)
          x=self.taCPMG
          print ('aqtime', x)
      if self.dataType=='Diffusion':
          if len(self.gradAmpArray) != nspectra:
              self.message('Data does not match gradient amplitude array')
             # return
          #self.bValueArray=(self.gradAmpArray*self.bcalfactor)**2
          x=self.bValueArray[:nspectra]
      if self.dataType=='RFattn':   #Rf attenuation data where RF attenuation is stepped acording to table RFattnDB
          try:
              self.RFattn=self.tntfile.DELAY['RFattnDB']
              logP=(self.RFampgain-self.RFattn)/10      #log power in mW, 
              #x=10**logP[0:nspectra]/1000     #chop off some of the array since the number of spectra may be less thane the table length, conver to Watts
              x= self.RFattn[0:nspectra]     #chop off some of the array since the number of spectra may be less thane the table length
          except:
             raise 
             self.message('Cannot find RF attenuation table')
             x=np.arange(nspectra)

      for j in range (nrepeats): #integrate FID/spectra
        for i in range(nspectra):
            self.dataInt[i,j]=np.trapz(self.data[:,i,j,0].real[self.iStart:self.jStop])
            #self.message("spectra {} Integrated value={:.3e}".format(j, self.dataInt[i,j]))
            #self.message ("{:.4e}".format(x[i]) + "   " + "{:.4e}".format(self.dataInt[i,j]))
      self.dataInt=self.dataInt/np.max(self.dataInt)        #Normalize integrated data to first spectrum
      if self.dataType=='Diffusion' and self.doNotFitDBelowSNR:     #if SNR flag then remove data points where the signal is below the SNR limit 1/self.SNRminfit
          self.snrClip=np.where(self.dataInt < 1/self.SNRminfit)
          x=np.delete(x,self.snrClip )
          self.bValueArray=np.delete(self.bValueArray,self.snrClip )
          self.clippedDataInt=np.zeros((x.shape[0],nrepeats))
          self.message('Clipping points below SNR from {} points to {} points'.format(nspectra,x.shape[0]), color='orange')
          for j in range (nrepeats):
              self.clippedDataInt[:,j] = np.delete(self.dataInt[:,j], self.snrClip)
          self.dataInt=self.clippedDataInt
      self.pw=self.plotWindow(self)
      self.pw.win.setWindowTitle('Integrated ' + self.dataType + ' Data repeat= ' + str(self.repeatIndex)+ ', ' + self.fileName)
      self.pw.show()
      if self.dataType=='FID':
        self.pw.setLogMode(logx=False,logy=False)
        self.pw.dplot.plot(x,self.dataInt[:,self.repeatIndex], pen=None, symbol='o',symbolSize=14,)
        self.pw.dplot.setLabel('bottom', "FID Index", **self.labelStyle) 
      if self.dataType=='HP':
        self.pw.setLogMode(logx=False,logy=False)
        xon=x[self.HPMWmask==0]
        yon=self.dataInt[:,self.repeatIndex][self.HPMWmask==0]
        self.pw.dplot.plot(xon,yon, pen=None, symbol='o',symbolSize=14,symbolBrush=(100,250,100))
        xoff=x[self.HPMWmask==1]
        yoff=self.dataInt[:,self.repeatIndex][self.HPMWmask==1]
        self.pw.dplot.plot(xoff,yoff, pen=None, symbol='o',symbolSize=14,symbolBrush=(250,100,100))
        self.pw.dplot.setLabel('bottom', "FID Index", **self.labelStyle)
        self.HPindex=np.average(yon)/ np.average(yoff)
        self.message("MW on={:.3f}, MW off={:.3f}, HP Enhancement={:.3f}".format(np.average(yon),np.average(yoff), self.HPindex))
        self.pw.dplot.setTitle("MW on={:.3f}, MW off={:.3f}, HP Enhancement={:.3f}".format(np.average(yon),np.average(yoff), self.HPindex), size='24pt', color=(100,200,200))
      if self.dataType=='IR':
        self.pw.setLogMode(logx=True,logy=False)
        header='IR data \n' +self.fileName + '\n TI(s) Integrated Signal'        
        for i in range(nrepeats):
          self.pw.dplot.plot(x, self.dataInt[:,i], pen=None, symbol=self.symb[i],symbolSize=24,border=self.kpen)
          self.pw.addData(x, self.dataInt[:,i], header=header)
          self.pw.addSNRData(x, self.SNR[:,i]) 
        self.pw.dplot.setLabel('bottom', "TI",units='s') #, **self.labelStyle)
        
      if self.dataType=='T1rho':
        self.pw.setLogMode(logx=True,logy=False)
        header='T1rho data \n' +self.fileName + '\n Integrated Signal vs spinlock time'        
        for i in range(nrepeats):
          self.pw.dplot.plot(x, self.dataInt[:,i], pen=None, symbol=self.symb[i],symbolSize=24,)
          self.pw.addData(x, self.dataInt[:,i], header=header)
          self.pw.addSNRData(x, self.SNR[:,i]) 
        self.pw.dplot.setLabel('bottom', "Spin lock time",units='s') #, **self.labelStyle)    
      if self.dataType=='RFattn':
        self.pw.setLogMode(logx=False,logy=False)
        self.pw.dplot.plot(x,self.dataInt[:,self.repeatIndex], pen=None, symbol='o',symbolSize=14,)
        self.pw.dplot.setLabel('bottom', "RF attenuation (dB)", **self.labelStyle)         
      if self.dataType=='nutation':
        ts=float(self.ui.leNutationStart.text())*1E-6
        ti=float(self.ui.leNutationIncrement.text())*1E-6
        self.tn=np.arange(nspectra)*ti+ts
        self.pw.setLogMode(logx=False,logy=False)
        self.pw.dplot.plot(self.tn,self.dataInt[:,self.repeatIndex], pen=None, symbol='o',symbolSize=14,)
        self.pw.dplot.setLabel('bottom', "tRF",units='s', **self.labelStyle) 
      if self.dataType=='ECCringdown':
        self.pw.setLogMode(logx=True,logy=False)
        self.pw.dplot.plot(self.gradRingdownDelay,self.FreqArray[:,self.repeatIndex], pen=None, symbol='o',symbolSize=14,)
        self.pw.dplot.setLabel('bottom', "Ringdown Time",units='s', **self.labelStyle) 
      if self.dataType=='CPMG':
        self.pw.setLogMode(logx=False,logy=True)
        header='CPMG data \n' +self.fileName + '\n tauCP(s) Integrated Signal'
        for i in range(nrepeats):
          self.pw.dplot.plot(x, self.dataInt[:,i], pen=None, symbol=self.symb[i],symbolSize=14,)
          self.pw.addData(x, self.dataInt[:,i], header=header)
          self.pw.addSNRData(x, self.SNR[:,i])
        self.pw.dplot.setLabel('bottom', "TCPMG",units='s', **self.labelStyle) 
      if self.dataType=='CPMG_tau_array':
        self.pw.setLogMode(logx=False,logy=True)
        header='CPMG tau array data \n' +self.fileName + '\n tauCP(s) Integrated Signal'
        for i in range(nrepeats):      
          self.pw.dplot.plot(x[:,i], self.dataInt[:,i], pen=None, symbol=self.symb[i],symbolSize=14,)
          self.pw.addData(x[:,i], self.dataInt[:,i], header=header)
          #self.pw.addSNRData(x[:,i], self.SNR[:,i])
        self.pw.dplot.setLabel('bottom', "Ta CPMG (s)")#, **self.labelStyle)    
      if self.dataType=='SE':
        self.pw.setLogMode(logx=False,logy=True)
        header='SE data \n' +self.fileName + '\n  Integrated Signal'
        for i in range(nrepeats):
            self.pw.dplot.plot(x, self.dataInt[:,i], pen=None, symbol=self.symb[i],symbolSize=14,)
            self.pw.addData(x, self.dataInt[:,i], header=header)
            self.pw.addSNRData(x, self.SNR[:,i])
        self.pw.dplot.setLabel('bottom', "TSE",units='s', **self.labelStyle)   
      if self.dataType=='Diffusion':
        if self.PGSEdelta>0:
            self.pw.dplot.setLabel('bottom',"b-Value (s/mm^2)") #, **self.labelStyle)
            self.pw.dplot.setLabel('left',"Normalized Integrated Signal")
        else:
            self.pw.dplot.setLabel('bottom',"Grad^2")
            self.pw.dplot.setLabel('left',"Normalized Integrated Signal")
        self.pw.setLogMode(logx=False,logy=True)
        header='Diffusion data \n' +self.fileName + '\n b(s/mm^2) Integrated Signal'
        for i in range(nrepeats):
          p=pg.mkPen(pg.intColor(i+self.penstep), width=1)  
          self.pw.dplot.plot(x, self.dataInt[:,i], pen=None, symbolBrush=pg.intColor(i), symbol=self.symb[i],symbolSize=12,)
          self.pw.addData(x, self.dataInt[:,i], header=header)
          self.pw.addSNRData(x, self.SNR[:,i])
   
  def plotSignalMax(self):
      '''Plots maximum of the magnitude of the Signal values'''
      self.pw=self.plotWindow(self)
      self.pw.win.setWindowTitle('FID Max ' + self.dataType + ' Data repeat= ' + str(self.repeatIndex)+ ', ' + self.fileName)
      self.pw.dataType='SignalMax'
      self.pw.show()
      nspectra=int(self.ui.lenSpectra.text())
      nrepeats=self.data.shape[2]
      self.dataMax=np.zeros((nspectra,nrepeats))        #array of maximum FID values
      x=self.bValueArray[:nspectra]
      for j in range (nrepeats): #integrate FID/spectra
        for i in range(nspectra):
            self.dataMax[i,j]=np.max(np.absolute(self.data[:,i,j,0]))
      self.dataMax= self.dataMax/np.max(self.dataMax)
      if self.dataType=='Diffusion':
        if self.PGSEdelta>0:
            self.pw.dplot.setLabel('bottom',"b-Value (s/mm^2)") #, **self.labelStyle)
        else:
            self.pw.dplot.setLabel('bottom',"Grad^2", units='au', unitPrefix=None)
        self.pw.setLogMode(logx=False, logy=True)
        header='Diffusion data \n' +self.fileName + '\n b(s/mm^2) Integrated Signal'
        for i in range(nrepeats):
          p=pg.mkPen(pg.intColor(i+self.penstep), width=1)  
          self.pw.dplot.plot(x, self.dataMax[:,i], pen=None, symbolBrush=pg.intColor(i), symbol=self.symb[i],symbolSize=12,)
          self.pw.addData(x, self.dataMax[:,i], header=header)

          
  def FFTData(self):
      '''does a discrete fast fourier transfrom on self.data, shift so 0 frequency is in the center''' 
      self.FIDData=copy.copy(self.data)
      s=self.data.shape
      self.ftData=np.zeros((s[0],s[1],s[2],s[3]),dtype='complex64')
      for j in range(self.nRepeats):
        for i in range(self.nSpectra):
            self.ftData[:,i,j,0]=fft.fftshift(fft.fft(self.data[:,i,j,0]))
      self.data=self.ftData
      self.ui.cbDataType.setCurrentIndex(self.ui.cbDataType.findText('Spectra (Hz)'))
      self.plotData()

  def inputData(self):
          self.message('<b>t90</b>(&mu;s)=' + self.ui.let90.text() +', <b>t180</b>(&mu;s)=' + self.ui.leT180.text(), report=True)
      

  def resetData(self):
    '''resets data to the original file data'''
    self.data=self.tntfile.DATA
    self.tfid=self.tntfile.fid_times()      #time points for FID waveforms
    self.freq=self.tntfile.freq_Hz()        #frequency points for spectra
    self.nRepeats=self.data.shape[2]
    self.ui.lenTD3.setText(str(self.nRepeats))    
    self.ui.sbRepeat.setMaximum(self.nRepeats-1)
    self.ui.cbDataType.setCurrentIndex(self.ui.cbDataType.findText('FID (s)'))
    self.Phase=0
    self.ui.hsPhase.setValue(self.Phase*10)
    self.ui.lePhase.setText("{:.1f}".format(self.Phase))
    self.plotData()
    
  def invertData(self):
    self.data=-self.data
    self.plotData()
    
  def removeLeadingTrailingPoints(self):
    '''trims first and last points from NMR FID or spectra to remove digital fileter artifacts or intial transients'''
    np, OK= QInputDialog.getInt(self,"trim input data", "# of points to remove from leading/trailing waveform", value=20)
    if OK:
        self.data=self.data[np:-np,:,:]
        self.tfid=self.tfid[np:-np]
        self.freq=self.freq[np:-np]
        self.nPoints=self.data.shape[0]
        self.ui.lenPoints.setText(str(self.nPoints))
        self.plotData()

  def subtractReferenceCurve(self):
    rc, OK= QInputDialog.getInt(self,"Subtract reference curve", "enter reference curve", value=1)
    if OK:
      refCurve=self.data[:,rc,:,:]  
      for i in range(self.nSpectra):
        self.data[:,i,:,:]-=refCurve
    self.plotData()
                      
  def deleteCurrentRepeat(self):
    '''resets data to original file data, deletes a the selected repeat to eliminate bad data or to analyze just one set'''
    self.repeatIndex=int(self.ui.sbRepeat.value())
#    self.data=np.delete(self.tntfile.DATA, self.repeatIndex,axis=2)  #input 4 dimensional data array, usually 4th dimension not used
    self.data=np.delete(self.data, self.repeatIndex,axis=2)  #input 4 dimensional data array, usually 4th dimension not used
    self.nRepeats=self.data.shape[2]
    self.ui.lenTD3.setText(str(self.nRepeats))    
    self.ui.sbRepeat.setMaximum(self.nRepeats-1)
    self.ui.cbDataType.setCurrentIndex(self.ui.cbDataType.findText('FID (s)'))
    self.plotData()
  
  def enterTIList(self):
    text, ok = QInputDialog.getText(self, 'Input TI', 
            'Enter TI list', text = str(self.TI).strip('[]'))
    if ok:
          if 's' in str(text):
              seper='s'
          else:
               seper = ' '
          self.TI=np.fromstring(str(text), sep=seper)
          self.message('TI(s) = ' + str(self.TI))
          
  def enterTEList(self):
    text, ok = QInputDialog.getText(self, 'Input TE', 
            'Enter TE list', text = str(self.TI).strip('[]'))
    if ok:
          if 's' in str(text):
              seper='s'
          else:
               seper = ' '
          self.TE=np.fromstring(str(text), sep=seper)
          self.message('TE(s) = ' + str(self.TE)) 
                     
  def enterCPMGList(self):
    text, ok = QInputDialog.getText(self, 'Input CPMG', 
            'Enter rephasing loop table (will add 1)', text=str(self.cpmgLoopTable).strip('[]'))
    if ok:
        self.cpmgLoopTable=np.fromstring(str(text), sep=' ')
        self.cpmgLoopIt=self.cpmgLoopTable+1
        self.message('CPMG loop list='  + str(self.cpmgLoopIt))
        
#*****************Nutation****************            
  def fitNutation(self):
      self.fitNutationData(self.tn,self.dataInt[:,0])  

  def fitNutationData(self,tn,data):
      """Fits nutation data with sinusoid, calls fitting routines in dampedSin Signal = a*(sin(pi/2*t/t90) * np.exp(-t/tau))"""
      self.fitx =np.arange(self.nfitpoints) * np.amax(tn) * 1.1 /self.nfitpoints  #generate RF pulse times for fit plot
      self.fity=np.zeros(self.nfitpoints)
      self.message ('<b>damped Sin 3-param fit</b>: t90(micros)=time for 90 deg flip, A=signal amplitude, tau=damping time')
      params=dampedSin.initialize (t=tn, s=data)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      fitoutput = lmfit.minimize(dampedSin.dSin,pdicti,args=(tn,data))
      pdict=fitoutput.params
      self.fity= dampedSin.dSin(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.message(lmfit.fit_report(pdict)+'\n')   #add complete fitting report to output report string
      self.pw.dplot.plot(self.fitx,self.fity, pen=self.bpen)
      self.t90Cal=float(pdict['t90'])
      self.NutTau=float(pdict['tau'])
      self.pw.addFitData(self.fitx,self.fity,title='t90: ' + "{:.6f}".format(self.t90Cal)+ "{:.6f}".format(self.NutTau))
      self.pw.addResData(tn,fitoutput.residual)
      self.pw.plotTitle = 't90(&mu;s)=' + "{:.4f}".format(1E6*self.t90Cal)+ ',  tau(&mu;s)=' + "{:.4f}".format(1E6*self.NutTau) + '; '
      self.pw.dplot.setTitle(self.pw.plotTitle)
      self.message("t90(us)= {:.2f}, B1(uT)={:.2f}, tau(us)={:.2f}".format(1E6*self.t90Cal,1E6*np.pi/(2*self.Gamma*self.t90Cal),1E6*self.NutTau), color='green', bold=True)
      cbText="{:.2f} {:.2f} {:.2f}".format(1E6*self.t90Cal,1E6*np.pi/(2*self.Gamma*self.t90Cal),1E6*self.NutTau)
      clipboard.setText(cbText)
      self.message('text sent to clipboard: ' + cbText, color='blue')
      return(fitoutput)  

  def plotNutationData(self,tn,data):
      self.RFamplitude=float(self.ui.leRFamplitude.text())*1.0E-6
      self.coilSensitivity()        #recalculates the coil sensitivity function
      self.nutateSignal=np.zeros(self.nfitpoints)
      t=np.arange(self.nfitpoints)*tn.max()/self.nfitpoints
      for i,tau in enumerate(t):
          self.nutateSignal[i]=nutatefit.signal(self.b1z,self.RFamplitude, tau)
      s=self.nutateSignal.max()
      self.nutateSignal=self.nutateSignal*self.dataInt.max()/s
      self.pw.dplot.plot(t,self.nutateSignal, pen=self.rpen)
      
  def showDataWindow(self): 
      self.pw=self.plotWindow(self)
      self.pw.win.setWindowTitle('Data Window ') 
      self.pw.show()

#**************T1,T2 Fitting******************************************************************    
  def makeFakeT1IRData(self): 
      self.pw=self.plotWindow(self)
      self.pw.win.setWindowTitle('Fake T1IR data= ') 
      self.pw.show()
      nspectra=int(self.ui.lenSpectra.text())
      nrepeats=int(self.ui.lenTD3.text())
      self.nRepeats = nrepeats
      self.dataInt=np.zeros((nspectra,nrepeats))
      text, ok = QInputDialog.getText(self, 'Input T1 (ms)', 'Enter T1(ms) list e.g.: 40,5,1')
      if ok:
            pass
      T1a=0.042
      T1b=0.005
      fT1b=0.1
      noise=0.02
      self.pw.setLogMode(True, False)
      self.pw.dplot.setLabel('bottom', "TI",units='s', **self.labelStyle)
      header=' Synthesized IR data \n' + '\n TI(s) Integrated Signal'      
      for i in range(nrepeats):
          self.dataInt[:,i]=((1-2*np.exp(-self.TI/T1a)) +fT1b*(1-2*np.exp(-self.TI/T1b)) + np.random.normal(scale=noise, size=nspectra))/(1+fT1b)
          self.pw.dplot.plot(self.TI, self.dataInt[:,i], pen=None, symbol='o',symbolSize=24,)
          self.pw.addData(self.TI, self.dataInt[:,i], header=header)
           
  def fitT1IR(self):
      self.nRepeats=self.dataInt.shape[1]
      self.T1fit=np.zeros(self.nRepeats)
      self.T1fitStdErr=np.zeros(self.nRepeats)
      self.T1deltafit=np.zeros(self.nRepeats)
      self.expected=np.zeros((self.dataInt.shape[0],self.dataInt.shape[1]))
      self.residualSTD=0
      self.plotTitle="T1-IR: "      
      for i in range(self.nRepeats):    #fit each repeat separately
          self.t1FitResults=self.fitT1IRData(self.TI,self.dataInt[:,i], repeat=str(i))        
          self.T1fit[i]=self.t1FitResults.params['T1'].value
          self.T1fitStdErr[i]=self.t1FitResults.params['T1'].stderr
          self.T1deltafit[i]=self.t1FitResults.params['delta'].value
          pdict=self.t1FitResults.params
          self.residualSTD+=np.std(self.t1FitResults.residual)
      self.residualSTD=self.residualSTD/self.nRepeats
      self.T1reported=np.mean(self.T1fit)
      self.T1StdErrReported=np.mean(self.T1fitStdErr)
      self.deltaReported=np.mean(self.T1deltafit)
      self.message('<b>Ave of fits, T1 mean and sd (ms)=</b>'  + "{:.6f}".format(self.T1reported*1000.0) + ' , ' + "{:.6f}".format(np.std(self.T1fit)*1000.0) +
            ' <b>delta mean and sd=</b>'  + "{:.6f}".format(self.deltaReported) + ' , ' + "{:.6f}".format(np.std(self.T1deltafit)) + ', Standard deviation of residuals=' + "{:.4e}".format(self.residualSTD)) 
      #fit all data together
      t1=np.repeat(self.TI, self.dataInt.shape[1])
      data=self.dataInt.flatten()
      self.t1FitResults=self.fitT1IRData(t1,data, repeat='all')        
      self.T1fitAll=self.t1FitResults.params['T1'].value
      self.T1deltafitAll=self.t1FitResults.params['delta'].value
      pdict=self.t1FitResults.params
      self.residualSTD=np.std(self.t1FitResults.residual)
      for i in range(self.nRepeats):
          self.expected[:,i]= T1IR.T1IR(pdict, self.TI, np.zeros(len(self.TI)))
      self.message('<b>Fitting All, T1 mean (ms)=</b>'  + "{:.6f}".format(np.mean(self.T1fitAll)*1000.0)  + 
                     ', <b>delta mean=</b>'  + "{:.6f}".format(np.mean(self.T1deltafitAll)) )
      par=self.lackOfFitTest(self.dataInt,self.expected, 3) #lack of fit test returns F statistic and p value
      self.T1IRFstat=par[0]
      self.T1IRpvalue=par[1]
      if self.T1IRpvalue>self.lofcriteria or self.residualSTD < self.maxResidualSTD:
        fittest=self.formatText('Quality of fit OK', color='green', bold=True)
      else:
        fittest=self.formatText('***Warning quality of T1-IR fit does not pass p-value test ***', color='red', bold=True)
      fitStats='Fstat='  + "{:.4f}".format(self.T1IRFstat ) + ', p-value=' + "{:.3e}".format( self.T1IRpvalue) + ', Standard deviation of residuals=' + "{:.4e}".format(self.residualSTD) + ',  ' + fittest
      self.message(fitStats)
      self.pw.addFitStats(fitStats)
      
      if  self.recoveryTime > self.T1reported*self.nT1sRequired:
        recoverytimetest=self.formatText('Recovery time,{:.3f}s, greater than {:.1f} T1s, sequence OK'.format(self.recoveryTime,self.nT1sRequired), color='green', bold=True)
      else:
        recoverytimetest=self.formatText('***Warning recovery time {:.3f} not adequate ***'.format(self.recoveryTime), color='red', bold=True)
      self.message(recoverytimetest)
      self.pw.addFitStats(recoverytimetest)

      if  self.deltaReported > self.deltaRequired:
        deltatest=self.formatText('Inversion efficiency delta greater than ' + str(self.deltaRequired) + ', sequence OK', color='green', bold=True)
      else:
        deltatest=self.formatText('***Inversion efficiency delta less than {:.3f} ,not adequate ***'.format(self.deltaRequired), color='red', bold=True)
      self.message(deltatest)
      self.pw.addFitStats(deltatest)
      for ii in range(self.nRepeats):
          self.message('T1r{} {:.5f} {:.5f}'.format(ii,self.T1fit[ii], self.T1deltafit[ii]))
      self.message(self.formatText('Acq Date, Temp(C), first peak FWHM(Hz) ={}   {:.3f}   {:.3f}'.format(self.tntfile.finish_time.isoformat(), self.temperature, self.FirstPeakfwhm) , color='blue', bold=False))      
      self.message(self.formatText('Reported T1 value(ms), delta ={:.3f}   {:.3f}'.format(self.T1reported*1000, self.deltaReported) , color='green', bold=True))
      cbText='{}  {:.3f}   {:.3f}  {:.3f}   {:.3f}'.format(self.tntfile.finish_time.isoformat(), self.temperature, self.FirstPeakfwhm, self.T1reported*1000, self.deltaReported)
      clipboard.setText(cbText)
      self.message('text sent to clipboard: ' + cbText, color='blue')

                      
  def fitT1IRData(self,TI,data, repeat='0'):
      """Fits T1-IR data, calls fitting routines in T1IR model = Si*(1-(1+delta)* np.exp(-TI/T1))"""
      self.fitx =np.arange(self.nfitpoints) * np.amax(TI) * 1.1 /self.nfitpoints + 0.001  #generate TIs for fit plot
      self.fity=np.zeros(self.nfitpoints)
      self.T1results=np.zeros(3)  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.message ('<b>T1-IR 3-param fit repeat=</b>'  + repeat + ': T1(s)=longitudinal spin relaxation time, A=signal at infinity, -delta=initial signal/signal at infinity')
      params=T1IR.initializeT1IR (TI=TI, data=data)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      fitoutput = lmfit.minimize(T1IR.T1IR,pdicti,args=(TI,data))
      pdict=fitoutput.params
      self.fity= T1IR.T1IR(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.message(lmfit.fit_report(pdict)+'\n')   #add complete fitting report to output report string
      self.pw.dplot.plot(self.fitx,self.fity, pen=self.bpen)
      self.pw.addFitData(self.fitx,self.fity,title='T1IR: ' + "{:.6f}".format(float(pdict['T1']))+ "{:.6f}".format(float(pdict['delta'])))
      self.pw.addResData(TI,fitoutput.residual)
      self.pw.plotTitle +='#' +str(repeat) + ' T1(ms)=' + "{:.3f}".format(1000*float(pdict['T1']))+ ',  &delta;=' + "{:.3f}".format(float(pdict['delta'])) + '; '
      self.pw.setPlotTitle()
      self.pw.dplot.setLabel('bottom', "TI (s)") #**self.labelStyle)
      self.pw.dplot.setLabel('left', "Signal") #, **self.labelStyle)
      return(fitoutput)  
 
  def fitCPMGT2(self):
      ''' Fits CPMG and SE data sets to get T2 by looping over repeats and calling CPMGT2data'''
      self.T2fit=np.zeros(self.nRepeats)
      nspectra=int(self.ui.lenSpectra.text())
      nspectra=self.dataInt.shape[0]
      nrepeats=self.dataInt.shape[1]
      self.expected=np.zeros((nspectra,self.dataInt.shape[1]))
      self.residualSTD=0
      if self.dataType=='CPMG':
          self.plotTitle="T2-CPMG: "
          ta=self.taCPMG[:nspectra]
      if self.dataType=='SE':
          self.plotTitle="T2-SE: "
          ta=self.tSE[:nspectra]
      #fits each repeat data set independently    
      for i in range(self.nRepeats):
          self.t2FitResults=self.fitCPMGT2Data(ta,self.dataInt[:,i], repeat=str(i))        
          self.T2fit[i]=self.t2FitResults.params['T2'].value
          pdict=self.t2FitResults.params
          self.residualSTD+=np.std(self.t2FitResults.residual)
          self.expected[:,i]= T2CPMG.T2CPMG(pdict, ta, np.zeros(nspectra))
      self.residualSTD=self.residualSTD/self.nRepeats
      self.T2reported=np.mean(self.T2fit)
      self.message('T2values(ms)={}'.format(np.array2string(self.T2fit*1000.0, precision=2)))
      self.message('<b>Ave of fits</b>, T2 mean and sd (ms)='  + "{:.4f}".format(self.T2reported*1000.0) + ' , '+ "{:.4f}".format(np.std(self.T2fit)*1000.0))
    #  fit all curves simulaneously         
      t2=np.repeat(ta, self.nRepeats)
      data=self.dataInt.flatten()
      self.t2FitResults=self.fitCPMGT2Data(t2,data, repeat='all')        
      self.T2fit=self.t2FitResults.params['T2'].value
      pdict=self.t2FitResults.params
      self.residualSTD=np.std(self.t2FitResults.residual)
      for i in range(self.nRepeats):      #calculates expected values for statistics
          self.expected[:,i]= T2CPMG.T2CPMG(pdict, ta, np.zeros(len(ta)))
      self.message('<b>Fitting All</b>')
      par=self.lackOfFitTest(self.dataInt,self.expected, 2)
      self.cpmgFstat=par[0]
      self.cpmgpvalue=par[1]
      if self.cpmgpvalue>self.lofcriteria or self.residualSTD<self.maxResidualSTD:
        fittest=self.formatText('quality of fit OK', color='green', bold=True)
      else:
        fittest=self.formatText('***Warning quality of fit does not pass p-value test ***', color='red', bold=True)
      fitStats='Fstat='  + "{:.4f}".format(self.cpmgFstat ) + ', p-value=' + "{:.3e}".format(self.cpmgpvalue) + ', Standard deviation of residuals=' + "{:.4e}".format(self.residualSTD) + ',  '+ fittest
      self.message(fitStats)
      self.pw.addFitStats(fitStats) 
      self.inhomogeneousFWHM=(1.0/self.T2star-1.0/self.T2reported)/np.pi
      if  self.inhomogeneousFWHM < self.maxInhLinewidth:
        fwhmtest='Inhomogeneous peak width ' + "{:.3f}".format(self.inhomogeneousFWHM) + '(Hz), is less than max allowed ' + "{:.3f}".format(self.maxInhLinewidth) + '(Hz),'
        fwhmtest+=self.formatText( 'Shim OK', color='green', bold=True)
      else:
        fwhmtest='***Warning Inhomogeneous peak width ' + "{:.3f}".format(self.inhomogeneousFWHM) + '(Hz), is greater than max allowed '+ "{:.3f}".format(self.maxInhLinewidth) + '(Hz),'
        fwhmtest+=self.formatText( 'Need to reshim', color='red', bold=True)
      self.message(fwhmtest)
      self.message(self.formatText('Acq Date, Temp(C), first peak FWHM(Hz) ={}   {:.3f}   {:.3f}'.format(self.tntfile.finish_time.isoformat(), self.temperature, self.FirstPeakfwhm) , color='blue', bold=False))      
      self.message(self.formatText('Reported T2 value (ms) =' +"{:.3f}".format(self.T2reported*1000) , color='green', bold=True))
      clipboard.setText('{:.3f}'.format(self.T2reported*1000))
      self.message('Text sent to clipboard: ' + '{:.3f}'.format(self.T2reported*1000), color='blue')

  def fitCPMGT2Data(self,TA,data, repeat='0'):
      """Fits T2CPMG data, calls fitting routines in T2CPMG model = Si*np.exp(-TA/T2)+B
      Returns fit output dictionary"""
      self.fitx =np.arange(self.nfitpoints) * np.amax(TA) * 1.1 /self.nfitpoints  #generate TIs for fit plot
      self.fity=np.zeros(self.nfitpoints)
      self.T2results=np.zeros(3)  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.message ('<b>T2-CPMG fit repeat=</b>' + repeat + ': T2(s)\n')
      params=T2CPMG.initializeT2CPMG (TE=TA, data=data)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      fitoutput = lmfit.minimize(T2CPMG.T2CPMG,pdicti,args=(TA,data))
      pdict=fitoutput.params
      self.fity= T2CPMG.T2CPMG(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.message(lmfit.fit_report(pdict)+'\n')   #add complete fitting report to output report string
      self.pw.dplot.plot(self.fitx,self.fity, pen=self.bpen)
      self.pw.addFitData(self.fitx,self.fity,title='T2CCPMG: ' + str(pdict['T2']))
      self.pw.addResData(TA,fitoutput.residual)
      self.pw.plotTitle +='R' +str(repeat) + ' T2(ms)=' + "{:.4f}".format(1000*float(pdict['T2'])) + '; '
      self.pw.setPlotTitle()
      self.pw.setWindowTitle('T2CPMG: ' + str(pdict['T2']))
      return(fitoutput)  

  def fitCPMGT2biEx(self):
      ''' Fits CPMG and SE data sets to biexponential model average over repeats and calling CPMGT2data'''
      self.T2afit=np.zeros(self.nRepeats)
      self.T2bfit=np.zeros(self.nRepeats)
      self.T2afitStdErr=np.zeros(self.nRepeats)
      self.T2bfitStdErr=np.zeros(self.nRepeats)
      self.Cafit=np.zeros(self.nRepeats)
      self.Cbfit=np.zeros(self.nRepeats)
      nspectra=int(self.ui.lenSpectra.text())
      nspectra=self.dataInt.shape[0]
      nrepeats=self.dataInt.shape[1]
      self.expected=np.zeros((nspectra,self.dataInt.shape[1]))
      self.residualSTD=0
      if self.dataType=='CPMG':
        xarray=self.taCPMG[0:nspectra]
      if self.dataType=='SE':
        xarray=self.tSE[0:nspectra]
      self.plotTitle="T2 bi-Exponential: "
      #fits each repeat data set independently    
      for i in range(self.nRepeats):
          self.t2FitResults=self.fitCPMGT2biExData(xarray,self.dataInt[:,i], repeat=str(i))        
          self.T2afit[i]=self.t2FitResults.params['T2a'].value
          self.T2bfit[i]=self.t2FitResults.params['T2b'].value
          self.T2afitStdErr[i]=self.t2FitResults.params['T2a'].stderr
          self.T2bfitStdErr[i]=self.t2FitResults.params['T2b'].stderr
          self.Cafit[i]=self.BiExpCa
          self.Cbfit[i]=self.BiExpCb
          pdict=self.t2FitResults.params
          self.residualSTD+=np.std(self.t2FitResults.residual)
          self.expected[:,i]= T2CPMGbiExp.ObjFnc(pdict, xarray, np.zeros(nspectra))
      self.residualSTD=self.residualSTD/self.nRepeats
      if np.mean(self.Cafit)>np.mean(self.Cbfit):   #order T2s by largest to smallest fraction
          self.T2areported=np.mean(self.T2afit)
          self.T2breported=np.mean(self.T2bfit)
          self.T2aStdErrReported=np.mean(self.T2afitStdErr)
          self.T2bStdErrReported=np.mean(self.T2bfitStdErr)
          self.CaReported=np.mean(self.Cafit)
          self.CbReported=np.mean(self.Cbfit)
      else:
          self.T2areported=np.mean(self.T2bfit)
          self.T2breported=np.mean(self.T2afit)
          self.T2aStdErrReported=np.mean(self.T2bfitStdErr)
          self.T2bStdErrReported=np.mean(self.T2afitStdErr)
          self.CaReported=np.mean(self.Cbfit)
          self.CbReported=np.mean(self.Cafit)
      self.message('<b>Ave of fits</b>, T2major mean and sd (ms)='  + "{:.4f}".format(self.T2areported*1000.0) + ' , '+ "{:.4f}".format(np.std(self.T2afit)*1000.0))
      self.message('<b>Ave of fits</b>, T2minor mean and sd (ms)='  + "{:.4f}".format(self.T2breported*1000.0) + ' , '+ "{:.4f}".format(np.std(self.T2bfit)*1000.0))
#      self.message('<b>T2short fraction, T2long fraction</b>, ='  + "{:.4f}, {:.4f}".format(self.CaReported, self.CbReported))
      self.message("<b>T2a(ms) Ca(%) T2b(ms) Cb(%)=</b>{:.4f} {:.4f} {:.4f} {:.4f}".format(self.T2areported*1000.0,self.CaReported*100,self.T2breported*1000.0, self.CbReported*100)) 

  def fitCPMGT2biExData(self,TA,data, repeat='0'):
      """Fits T2CPMG and T2SE data, calls fitting routines in T2CPMG bi exponential model = Si*np.exp(-TA/T2)+Sib*np.exp(-TA/T2b) +B"""
      self.fitx =np.arange(self.nfitpoints) * np.amax(TA) * 1.1 /self.nfitpoints  #generate TIs for fit plot
      self.fity=np.zeros(self.nfitpoints)
      self.T2results=np.zeros(3)  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.message ('<b>T2-CPMG biExp fit repeat=</b>' + repeat + ': T2(s)\n')
      params=T2CPMGbiExp.initialize (TE=TA, data=data)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      fitoutput = lmfit.minimize(T2CPMGbiExp.ObjFnc,pdicti,args=(TA,data))
      pdict=fitoutput.params
      self.fity= T2CPMGbiExp.ObjFnc(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.message(lmfit.fit_report(pdict)+'\n')   #add complete fitting report to output report string
      self.pw.dplot.plot(self.fitx,self.fity, pen=self.bpen)
      self.BiExpSignalRatio=float(pdict['Sib'])/float(pdict['Si'])
      self.BiExpCa=float(pdict['Si'])/(float(pdict['Si'])+float(pdict['Sib']))
      self.BiExpCb=float(pdict['Sib'])/(float(pdict['Si'])+float(pdict['Sib']))
      self.pw.addFitData(self.fitx,self.fity,title='T2CCPMG biExp: ' + str(pdict['T2a']))
      self.pw.addResData(TA,fitoutput.residual)
      self.pw.plotTitle +='R' +str(repeat) + ' T2a(ms)=' + "{:.4f}".format(1000*float(pdict['T2a'])) + ' T2b(ms)=' + "{:.4f}".format(1000*float(pdict['T2b']))+ '; '+ ' Si/Sib=' + "{:.4f}".format(self.BiExpSignalRatio)
      self.pw.setPlotTitle()
      self.pw.setWindowTitle('T2 bi Exp: ' + str(self.fileName))
      return(fitoutput)
  
  def fitCPMGT2_tau_array(self):
      ''' Fits CPMG tau array data, outputs T2 vs refocus time'''
      self.T2fit=np.zeros(self.nRepeats)
      nspectra=int(self.ui.lenSpectra.text())
      nspectra=self.dataInt.shape[0]
      nrepeats=self.dataInt.shape[1]
      self.expected=np.zeros((nspectra,self.dataInt.shape[1]))
      self.residualSTD=0
      self.plotTitle="T2-CPMG tau array: "
      ta=self.taCPMG
      #fits each repeat data set independently
      summary=''    
      for i in range(self.nRepeats):
          self.t2FitResults=self.fitCPMGT2Data(ta[:,i],self.dataInt[:,i], repeat=str(i))        
          self.T2fit[i]=self.t2FitResults.params['T2'].value
          pdict=self.t2FitResults.params
          self.residualSTD+=np.std(self.t2FitResults.residual)
          self.expected[:,i]= T2CPMG.T2CPMG(pdict, ta[:,i], np.zeros(nspectra))
 
      self.residualSTD=self.residualSTD/self.nRepeats
      self.T2reported=np.mean(self.T2fit)
      self.message('<b>T2 CPMG tau array</b>')
      self.message('Tau(ms) RefocusFrequency(Hz), T2(ms)')
      for  i in range(self.nRepeats):        
          self.message("{:.2f}  {:.1f}  {:.1f}".format((self.cpmg_tau_array[i]+self.tauCPMG)*1000.0, 0.5/(self.cpmg_tau_array[i]+self.tauCPMG), self.T2fit[i]*1000.0))
   
      
  def fitT1rho(self):
      ''' Fits T1Rho data with simple exponential'''
      self.T1rhofit=np.zeros(self.nRepeats)
      nspectra=int(self.ui.lenSpectra.text())
      nspectra=self.dataInt.shape[0]
      nrepeats=self.dataInt.shape[1]
      self.expected=np.zeros((nspectra,self.dataInt.shape[1]))
      self.residualSTD=0
      xarray=self.TI[0:nspectra]
      self.plotTitle="T1 rho: "
      #fits each repeat data set independently    
      for i in range(self.nRepeats):
          self.T1rhoFitResults=self.fitT1rhoData(xarray,self.dataInt[:,i], repeat=str(i))        
          self.T1rhofit[i]=self.T1rhoFitResults.params['tau'].value
          pdict=self.T1rhoFitResults.params
          self.residualSTD+=np.std(self.T1rhoFitResults.residual)
          self.expected[:,i]= expFit.ObjFnc(pdict, xarray, np.zeros(nspectra))
      self.residualSTD=self.residualSTD/self.nRepeats
      self.T1rhoReported=np.mean(self.T1rhofit)
      self.message('<b>Ave of fits</b>, T1rho (ms)='  + "{:.4f}".format(self.T1rhoReported*1000.0) )


  def fitT1rhoData(self,t,data, repeat='0'):
      """Fits t1rho with exponetial"""
      self.fitx =np.arange(self.nfitpoints) * np.amax(t) * 1.1 /self.nfitpoints  #generate TIs for fit plot
      self.fity=np.zeros(self.nfitpoints)
      self.T1rhoresults=np.zeros(3)  #T1 fitting results, first index = ROI, second index = parameter referred to in T1Params
      self.message ('<b>T1 rho Exp fit repeat=</b>' + repeat + ': T1rho(s)\n')
      params=expFit.initialize (t=t, data=data)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      fitoutput = lmfit.minimize(expFit.ObjFnc,pdicti,args=(t,data))
      pdict=fitoutput.params
      self.fity= expFit.ObjFnc(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.message(lmfit.fit_report(pdict)+'\n')   #add complete fitting report to output report string
      self.pw.dplot.plot(self.fitx,self.fity, pen=self.bpen)
      self.pw.addFitData(self.fitx,self.fity,title='T1rho Exp: ' + str(pdict['tau']))
      self.pw.addResData(t,fitoutput.residual)
      self.pw.plotTitle +='R' +str(repeat) + ' T1rho(ms)=' + "{:.4f}".format(1000*float(pdict['tau'])) 
      self.pw.setPlotTitle()
      self.pw.setWindowTitle('T1rho Exp: ' + str(self.fileName))
      return(fitoutput)  
#*******************Diffusion calculations/fitting*************************************************

  def STbvalue(self,  g=0.1, delta=0.01, Delta=0.01, risetime=0.0001, pulsetype='trap'):
    '''calculates  diffusion b-values according to generalized Stejskal-Tanner equation
     delta=grad pulse duration, Delta=grad pulse separation, risetime is the the grad pulse rise and falltime, all times in s'''
    gamma=self.Gamma
    if pulsetype=='trap':
        epsilon=risetime/delta
        sigma=1-epsilon
        lamda=0.5
        kappa=0.5-sigma/6+epsilon**3/60/sigma**2-epsilon**2/12/sigma
    if pulsetype=='hSin':
        sigma=2/np.pi
        lamda=0.5
        kappa=3.0/8.0    
    b=1E-6*(sigma*gamma*delta*g)**2*(Delta-2*(lamda-kappa)*delta)       #value calculated in SI s/m^2; convert to s/mm^2, exception to SI rule!
    self.message('<b>b-Value calculated using ST formula:</b> {} pulse, risetime(ms)={:.4f}, duration(ms)={:.4f}, pulse spacing(ms){:.4f}'.format(pulsetype, risetime*1000, delta*1000, Delta*1000))
#    bRect=1E-6*(gamma*delta*g)**2*(Delta-delta/3)
    return b    
                
  def fitDiffusion(self, dataType='', dataStart=0, dataStop=-1, Kurtosis=False, biExp=False):
      ''' Fits Diffusion data sets to get ADC values by looping over repeats and calling fitDiffusiondata(),
      Will start at data point dataStart and end at datapoint dataStop
      Can use a simple exponential model or an exponential model with a kurtosis term'''
      self.ADCfit=np.zeros(self.nRepeats)   #largest diffusion component
      self.ADC2fit=np.zeros(self.nRepeats)  #minor diffusion component
      self.ADCminorFracfit=np.zeros(self.nRepeats) #Fraction of minor diffusion component
      self.Kfit=np.zeros(self.nRepeats)
      self.bMax=np.zeros(self.nRepeats) #maximum bvalue used in kurtosis fits
      self.bMaxk=np.zeros(self.nRepeats)    #maximum allowed bvalue for Kurtosis bMaxk=1/3KD
      nspectra=int(self.ui.lenSpectra.text())
      if dataType=='SignalMax':
          Ddata=self.dataMax[dataStart:dataStop,:]
      else:
          Ddata=self.dataInt[dataStart:dataStop,:]
      nspectra=Ddata.shape[0]
      nrepeats=Ddata.shape[1]
      self.expected=np.zeros((nspectra,Ddata.shape[1]))      #these are the expected values using the fit parameters
      self.residualSTD=0
      if Kurtosis:
          self.plotTitle="Diffusion Kurtosis-PGSE: "
      elif biExp:
          self.plotTitle="Diffusion Biexponetial-PGSE: "
      else:
          self.plotTitle="Diffusion Simple Exponential-PGSE: "
      #fits each repeat data set independently    
      for i in range(self.nRepeats):
          fit=self.fitDiffusionData(self.bValueArray[dataStart:dataStop],Ddata[:,i], repeat=str(i), Kurtosis=Kurtosis, biExp=biExp)
          self.DiffusionFitResults=fit[0]    
          self.ADCfit[i]=self.DiffusionFitResults.params['ADC'].value
          self.ADC2fit[i]=self.DiffusionFitResults.params['ADC2'].value
          print(self.ADC2fit[i])
          self.Kfit[i]=self.DiffusionFitResults.params['K'].value
          self.ADCminorFracfit[i]=self.DiffusionFitResults.params['Si2'].value/(self.DiffusionFitResults.params['Si'].value+self.DiffusionFitResults.params['Si2'].value)
          pdict=self.DiffusionFitResults.params
          self.residualSTD+=np.std(self.DiffusionFitResults.residual)
          self.expected[:,i]= DiffusionPGSE.model(pdict, self.bValueArray[dataStart:dataStop], np.zeros(len(self.bValueArray[dataStart:dataStop])))
          self.bMax[i]=fit[1]
          self.bMaxk[i]=fit[2]
      self.residualSTD=self.residualSTD/self.nRepeats
      self.ADCreported=np.mean(self.ADCfit)
      self.ADC2reported=np.mean(self.ADC2fit)
      self.Kreported=np.mean(self.Kfit)
      self.message('<b>Ave of fits</b>, ADC mean and sd (10^-3 mm^2/s)= {:.4f}, {:.4f}'.format(self.ADCreported*1000.0, np.std(self.ADCfit)*1000.0)) 
      if Kurtosis:
          self.message('<b>Ave of fits</b>, Kurtosis mean and sd= {:.4f}, "{:.4f}'.format(self.Kreported, np.std(self.Kfit)))
      if biExp:
          self.message('<b>Ave of fits</b>, minor ADC mean and sd= {:.4f}, "{:.4f}'.format(self.ADC2reported, np.std(self.ADC2fit)))
      #  fit all curves simulaneously      
      bValues=np.repeat(self.bValueArray[dataStart:dataStop], self.nRepeats)
      data=Ddata.flatten()
      fit=self.fitDiffusionData(bValues,data, repeat='all', Kurtosis=Kurtosis, biExp=biExp)
      self.DiffusionFitResults=fit[0]        
      self.ADCfitall=self.DiffusionFitResults.params['ADC'].value
      pdict=self.DiffusionFitResults.params
      self.residualSTD=np.std(self.DiffusionFitResults.residual)
      for i in range(self.nRepeats):      #calculates expected values for statistics
          self.expected[:,i]= DiffusionPGSE.model(pdict, self.bValueArray[dataStart:dataStop], np.zeros(len(self.bValueArray[dataStart:dataStop])))
      self.message('<b>Fitting All, ADCall (10^-3 mm^2/s)=</b>' +"{:.3f}".format(self.ADCfitall*1000))
      par=self.lackOfFitTest(Ddata,self.expected, 2)
      self.PGSEFstat=par[0]
      self.PGSEpvalue=par[1]
      if self.PGSEpvalue>self.lofcriteria or self.residualSTD<self.maxResidualSTD:
        fittest=self.formatText('quality of fit OK', color='green', bold=True)
      else:
        fittest=self.formatText('***Warning quality of fit does not pass p-value test ***', color='red', bold=True)
      fitStats='Fstat='  + "{:.4f}".format(self.PGSEFstat ) + ', p-value=' + "{:.3e}".format(self.PGSEpvalue) + ', Standard deviation of residuals=' + "{:.3e}".format(self.residualSTD) + ',  '+ fittest
      self.message(fitStats, color='blue')
      self.pw.addFitStats(fitStats) 
      self.message('Measured ADC values (10^-3 mm^2/s)=' + np.array2string(self.ADCfit*1000, precision=3))

      self.message(self.formatText('Reported diffusion coefficient (10^-3 mm^2/s) =' +"{:.3f}".format(self.ADCreported*1000) , color='green', bold=True))
      self.Dcorrected=self.ADCreported-self.eddyGxCorrection
      self.message(self.formatText('Eddy corrected {}-diffusion coefficient(10^-3 mm^2/s) = {:.3f}, ECC={:.3f}'.format(self.gradOrientation, self.Dcorrected*1000 ,self.eddyGxCorrection*1000), color='green', bold=True))
      if self.gradOrientation=='X':
          self.DxReported=self.Dcorrected
          self.Dx2Reported=self.ADC2reported
          self.KxReported=self.Kreported
      if self.gradOrientation=='Y':
          self.DyReported=self.Dcorrected
          self.Dy2Reported=self.ADC2reported
          self.KyReported=self.Kreported    
      if self.gradOrientation=='Z':
          self.DzReported=self.Dcorrected
          self.Dz2Reported=self.ADC2reported
          self.KzReported=self.Kreported
      if Kurtosis:
          self.message(self.formatText('Measured Kurtosis, bMax(s/mm^2), bMaxk(s/mm^2) =' + np.array2string(self.Kfit, precision=4)+ np.array2string(self.bMax, precision=4)+ np.array2string(self.bMaxk, precision=4), color='blue', bold=True))     
          
  def fitDiffusionData(self,bValue,data, repeat='0', Kurtosis=False, biExp=False):
      """Fits DiffusionPGSE data, calls fitting routines in DiffusionPGSE model = model = Si*np.exp(-b*D+K*x**2/6)+ Si2*np.exp(-b*D2)+B"""
      self.fitx =np.arange(self.nfitpoints) * np.amax(bValue) * 1.1 /self.nfitpoints  #generate bvalues for fit plot
      self.fity=np.zeros(self.nfitpoints)
      self.bValueresults=np.zeros(3)  #b-value fitting results
      if Kurtosis:
          self.message ('<b>Diffusion/Kurtosis-PGSE fit repeat=</b>' + repeat + ': bValue(s)\n')
      if biExp:
          self.message ('<b>Diffusion BiExponential-PGSE fit repeat=</b>' + repeat + ': bValue(s)\n')
      else:
        self.message ('<b>Diffusion-PGSE fit repeat=</b>' + repeat + ': bValue(s)\n')
      params=DiffusionPGSE.initialize (bvalue=bValue, data=data, Kurtosis=Kurtosis, biExp=biExp)
      pdicti=params[0] #parameter dictionary
      plist=params[1] #parameter list
      fitoutput = lmfit.minimize(DiffusionPGSE.model,pdicti,args=(bValue,data))
      pdict=fitoutput.params
      bmax=np.amax(bValue)
      bmaxk=0.0     #range overwhich kurtosis can be fit
      if Kurtosis:  #need to restrict the range of the fit
          bmaxk=3/(pdict['ADC'].value*pdict['K'].value) #do not fit for bvalues above this value
          ifit=self.fitx.shape[0]   #nummber points in the fit to show
          while bmax>bmaxk:

            for i,b in enumerate(bValue):
              if b>bmaxk:
                  bValue=np.delete(bValue,[i])
                  data=np.delete(data,[i])
            fitoutput = lmfit.minimize(DiffusionPGSE.model,pdicti,args=(bValue,data))
            pdict=fitoutput.params
            bmaxk=3/(pdict['ADC'].value*pdict['K'].value) #do not fit for bvalues above this value
            bmax=np.amax(bValue)
            ifit=np.argmax(self.fitx>bmaxk)
      self.fity= DiffusionPGSE.model(pdict, self.fitx, np.zeros(len(self.fitx)))
      self.message(lmfit.fit_report(pdict)+'\n')   #add complete fitting report to output report string
      self.pw.plotTitle +='R' +str(repeat) + ' ADC(10^-3 mm^2/s)=' + "{:.4f}".format(1000*float(pdict['ADC'])) + '; '         
      if Kurtosis:
          self.pw.dplot.plot(self.fitx[:ifit+1],self.fity[:ifit+1], pen=self.bpen)
          self.pw.plotTitle += '; Kurtosis={:.4f}, bmax={:.2f}, bmaxk={:.2f}'.format(pdict['K'].value,bmax,bmaxk) 
      else:
          self.pw.dplot.plot(self.fitx,self.fity, pen=self.bpen)
      self.pw.addFitData(self.fitx,self.fity,title='Diffusion PGSE: ADC= ' + str(pdict['ADC']))
      self.pw.addResData(bValue,fitoutput.residual)   
      self.pw.setPlotTitle()
      self.pw.setWindowTitle('Diffusion PGSE: ' + self.fileName)
      return([fitoutput, bmax, bmaxk]) 
  
  def diffusionSummary(self, kurtosis=False, biExp=False):
      '''Outputs summary of diffusion analysis assume x and y diffusion data have been taken, need to be generalized to 1, 2, or 3 axes'''
      self.message(self.formatText('Diffusion Summary ', bold=True, color='blue'))
      self.DavReported=(self.DxReported+self.DyReported)/2
      if self.T1reported==0:
        recoverytimetest=self.formatText('***Recovery time(s)= {:.3f} ?adequate? ***'.format(self.recoveryTime), color='orange', bold=True)
      elif self.recoveryTime>5*self.T1reported:
        recoverytimetest=self.formatText('***Recovery time(s)= {:.3f} adequate ***'.format(self.recoveryTime), color='green', bold=True)
      else:
        recoverytimetest=self.formatText('***Recovery time(s)= {:.3f} not adequate ***'.format(self.recoveryTime), color='red', bold=True)
      self.message(recoverytimetest)
      self.message('ImaxGx, ImaxGy(A)={:.3f} {:.3f}, Gradient Cal: Gxcal, Gycal(mT/m/A)={:.3f} {:.3f}'.format(self.Axmax,self.Aymax,1000*self.GxCal, 1000*self.GyCal), color='green', bold=True)
      self.message('Eddy Current Corrections (ECC): GxECC,GyECC(10^-3 mm^2/s)={:.4f} {:.4f}'.format(1000*self.eddyGxCorrection,1000*self.eddyGyCorrection), color='green', bold=True)
      self.message('Corrected Diffusion Coefficients: Dx, Dy, Dav(10^-3 mm^2/s)= {:.3f}  {:.3f} {:.3f} '.format(1000*self.DxReported, 1000*self.DyReported,1000*self.DavReported), color='green', bold=True)
      self.message('Diffusion Summary: FWHM(Hz), Recovery time, ImaxGx, ImaxGy, GxCal, GyCal, ECCx, ECCy, Dx, Dy, Dav') 
      self.message('{:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.4f} {:.4f} {:.3f}  {:.3f} {:.3f} '.format(self.FirstPeakfwhm, self.recoveryTime, self.Axmax,self.Aymax ,1000*self.GxCal, 1000*self.GyCal,1000*self.eddyGxCorrection,1000*self.eddyGyCorrection,1000*self.DxReported, 1000*self.DyReported,1000*self.DavReported), color='green', bold=True)
      if kurtosis or self.Kreported>0:
          self.Kav=(self.KxReported + self.KyReported)/2
          self.message(' x Kurtosis, y Kurtosis, average Kurtosis= {:.4f}, {:.4f}, {:.4f}'.format(self.KxReported, self.KyReported,self.Kav ), color='blue', bold=True)
      if biExp:
          biexmessage='  {:.4f}  {:.4f}  {:.4f}'.format(1000*self.Dx2Reported, 1000*self.Dy2Reported,1000*(self.Dx2Reported+self.Dy2Reported)/2)
          self.message( 'ADC2x, ADC2y, ADC2av='+ biexmessage, color='blue', bold=True)
      else:
          biexmessage=''
      cbText='{}  {:.3f}  {:.3f}   {:.3f}  {:.3f}   {:.3f}  {:.3f}   {:.4f}  {:.4f}   {:.3f}   {:.3f}   {:.3f}'\
                        .format(self.tntfile.finish_time.isoformat(), self.FirstPeakfwhm, self.recoveryTime, self.Axmax,self.Aymax ,1000*self.GxCal, 1000*self.GyCal,1000*self.eddyGxCorrection,\
                        1000*self.eddyGyCorrection,1000*self.DxReported, 1000*self.DyReported,1000*self.DavReported) + biexmessage
      clipboard.setText(cbText) 
      self.message('Text sent to clipboard: ' + cbText, color='blue')

  def makeFakeDiffusionData(self): 
      self.pw=self.plotWindow(self)
      self.pw.win.setWindowTitle('Fake Diffusion data= ') 
      self.pw.show()
      nspectra=int(self.ui.lenSpectra.text())
      nrepeats=int(self.ui.lenTD3.text())
      self.nRepeats = nrepeats
      self.dataInt=np.zeros((nspectra,nrepeats))
      Da=1E-3
      Db=0.5E-3
      fDb=0.2
      noise=0.002
      bmax=3000
      text, ok = QInputDialog.getText(self, 'Input D(mm^2/s)', 'Enter bmax, Da, Db, fb, 1/SNR',text=self.FakeDiffusionParameters)
      if ok:
            inp=np.fromstring(text, dtype=float, sep=',')
            self.FakeDiffusionParameters=text
      try:
          bmax=inp[0]
          Da=inp[1]
          Db=inp[2]
          fDb=inp[3]
          noise=inp[4]
      except:
          pass
      self.bValueArray=np.linspace(0,bmax, nspectra)
      self.pw.setLogMode(False,True)
      self.pw.dplot.setLabel('bottom', "b", **self.labelStyle)
      header=' Synthesized Diffusion data \n' + '\n b-value Integrated Signal'      
      for i in range(nrepeats):
          self.dataInt[:,i]=(1-fDb)*np.exp(-self.bValueArray*Da) + fDb*np.exp(-self.bValueArray*Db) + np.random.normal(scale=noise, size=nspectra)
          self.pw.dplot.plot(self.bValueArray, self.dataInt[:,i], pen=None, symbol='o',symbolSize=18,)
          self.pw.addData(self.bValueArray, self.dataInt[:,i], header=header)
      self.pw.fitDiffusion() 
      self.pw.win.setWindowTitle('Fake Diffusion data: bmax, Da, Db, fb, 1/SNR=' + text)
      self.message("Fake Diffusion data created: bmax, Da, Db, fb, 1/SNR= " + self.FakeDiffusionParameters, color='blue')            
      
  def useMeasuredBvalues(self):  #toggles between measured and calculated b-values
      self.buseMeasureBvalues=True
      self.bValueArray=self.measuredBValues
      self.message( '<b>Using measured b-values</b> ')
  
  def useCalculatedBvalues(self):  #toggles between measured and calculated b-values
      self.buseMeasureBvalues=False
      self.message( '<b>Using calculated b-values</b> ')
      self.resetGradandBvalues()        
  
  def resetGradandBvalues(self):
      '''read in gradinet parameters and recalculate gradient strength and bvalues'''
      self.message('<b>Recalculating gradient and b-values</b> ') 
      self.PGSEgrad=float(self.ui.leGradPW.text())/1000.0
      self.PGSEdelta=float(self.ui.lePGSEdelta.text())/1000.0
      self.Amax=float(self.ui.leImax.text())
      self.GCal=float(self.ui.leGradCal.text())/1000 
      #self.gradAmpArray=self.tntfile.DELAY[self.pgseLoopID].astype(float)     #tnt arrays can be interpreted as floats
#      self.message('<b>PGSE grad amplitude table:</b> ' + str(self.gradAmpArray))
      self.gradientCurrent=self.GrAmp/100*self.gradAmpArray/100*self.Amax
#      self.message('<b>PGSE Current Pulse(A):</b> ' + str(self.gradientCurrent))
      self.gradientAmplitude=self.gradientCurrent*self.GCal   #in T/m, Amax is current when DAC is set to 100
      self.message('<b>gradient strength (mT/m)=</b>' + str(self.gradientAmplitude*1000))
      if self.PGSEdelta>0 and self.PGSEgrad>0:
        gPulseWidth=self.PGSEgrad+2*self.gradPulseRiseTime
        gPulseDelay=gPulseWidth+2*self.PGSEdelta+2*self.t90+2*self.PGSEintergradTimeDelay
        self.bValueArray=self.STbvalue(g=self.gradientAmplitude, delta=gPulseWidth, Delta=gPulseDelay, risetime=self.gradPulseRiseTime, pulsetype=self.GradPulseType)
        #self.bValueArray=1E-6*(self.Gamma*gradAmp*self.PGSEgrad)**2*(2.0*self.PGSEgrad/3+2*self.PGSEdelta+2*self.t90+2*self.PGSEintergradTimeDelay)
      else:
        self.bValueArray=self.gradAmpArray**2
                #b-value in s/mm^2   so ADC in mm^2/s!!!!
      self.message( '<b>b-values (s/mm^2)=</b> ' + np.array2string(self.bValueArray, precision=2, separator=','))         #"{:.4f}".format(

  def updateGcalImax(self):
      '''Updates gradient calibration and max current from current settings in the text boxes'''
      self.Amax=float(self.ui.leImax.text())
      self.GCal=float(self.ui.leGradCal.text())/1000
      direction=self.ui.cbGradOrientation.currentText()
      if direction=='X':
            self.Axmax=self.Amax
            self.GxCal=self.GCal
      if direction=='Y':
            self.Aymax=self.Amax
            self.GyCal=self.GCal
      if direction=='Z':
            self.Azmax=self.Amax
            self.GzCal=self.GCal
      
      
  def inputGradientCurrentTraces(self):
      '''Inputs Gradient Waveform time,Gx,Gy,Gz; makes an array of peak currents used in calibrating gradients, calulates q and b-values'''
      #QMessageBox.about(self, "Input Gradient Waveforms", "1. Looks for gxcal, gycal, gzcal in file name to determine grad orientation,<br> 2. Assumes col0=time(s), col1=Gx(V), col2=Gy(V), col3=Gz(V),<br> 3. Requires first 10% to be zero to determine baseline, <br> 4. Finds 180 pulse time from filename key tflip=")
      self.message(self.formatText('Input Gradient Current Traces: ' + self.gradOrientation + '-direction', bold=True, color='blue'))
      f = QFileDialog.getOpenFileNames(self,"Open Gradient Current Traces " + self.gradOrientation+ '-direction',filter="Gradient Current Files (*.csv *.txt)")
      if type(f)==tuple:    #passes  Qstring with PyQt4 and a tuple of strings with PyQt5
        if f[1]=='':
          return None
        self.fileNames=f[0]
      else:
        if not f:
          return None
        self.fileNames=str(f)      #make sure fileName is a string, not a Qstring
      fn=str(self.fileNames[0]).lower()
      try:   #try to find the maximum DAC output in the filename, dacmax=max table gr0:2*GrAmp/100, info is not in scope files       
          dacmax=self.DACmaxString.findall(fn)[0].replace('dacmax=','')
          self.DACmax=float(dacmax)    #maximum gradient DAC output
      except:
          self.DACmax=0
      try:   #try to find the 180 pulse time in the filename       
          tflip=float(self.tflipString.findall(fn)[0].replace('tflip=',''))    #time in ms when 180 degree pulse is applied   
      except:
          tflip=0
      nFiles= len(self.fileNames)
      nwf=1 #waveform to analize
      if fn.find('_gx_') >=0 or fn.find('gxcal') >=0 or fn.find('gxmeas')>=0 or fn.find('gxdif')>=0 :      #Set gradient orientation  and waveform to be analyzed
            self.setGradOrientation('X')
            nwf=1
      if fn.find('_gy_') >=0 or fn.find('gycal') >=0 or fn.find('gymeas')>=0 or fn.find('gydif')>=0:
            self.setGradOrientation('Y')
            nwf=2
      if fn.find('_gz_') >=0 or fn.find('gzcal') >=0 or fn.find('gzmeas')>=0 or fn.find('gzdif')>=0:
            self.setGradOrientation('Z')
            nwf=3           
      self.measuredBValues=np.zeros(nFiles)
      self.gradQTraces=[]
      self.gradCurrentTraces=[]
      self.gradGTraces=[]
      self.message('<b>Input current traces:</b> ' + self.fileNames[0])
      self.message("Gradient monitor(A/V)={:.2f}, Gcal(mT/m/A)={:.2f}".format(self.cmCal, 1000*self.GCal))
      self.message("I(A), ave Gradient(mT/m), q/2&pi;(1/mm), b(s/mm^2)")
      self.measuredI=np.zeros(nFiles)
      if tflip>0:
          self.message('Flipping q after '+ str(tflip) + 'ms')
      for i in range(nFiles):
        fileName = str(self.fileNames[i])
        gradwf=np.genfromtxt(fileName,delimiter=',',skip_header=3)
        gtime=gradwf[:,0]
        if gradwf.shape[1]<nwf:
            nwf=1
        V=gradwf[:,nwf]
        if np.mean(V)>=0: #calculate amplitude of positive pulses
          Vm=V.max()
          Vmav=np.mean(V[V > 0.97*Vm])
        else:   #calculate amplitude of negative pulses
          Vm=V.min()
          Vmav=np.mean(V[V < 0.97*Vm])
        Ioffset=int(0.08*V.shape[0])  #use the first 8% of data to determine baseline
        V0=np.mean(V[0:Ioffset])
        V-=V0       #subtract voltage offset at 0 gradient pulse
        Imav=(Vmav-V0)*self.cmCal     #calculate average current during gradient pulse
        self.measuredI[i]=Imav
        gradient=V*self.cmCal*self.GCal        #Change current monitor voltage into field gradient in T/m
        self.gradCurrentTraces.append(np.array((gtime,gradient)).T)
        Q=self.Gammaf*integrate.cumtrapz(gradient,gtime, initial=0)/1E6
        q=np.array([gtime,Q])
        if tflip>0:
            i180 = (np.abs(gtime - tflip)).argmin()
            q180=q[1,i180]
            q[1,i180:]-=2*q180      #flip sign of q after 180 pulse
        else:
            q180=-1
        self.gradQTraces.append(q)
        q2=q[1,:]*q[1,:]*4*np.pi**2     #calculate q^2
        bvalue=np.trapz(q2, x=q[0,:])/1000      #integrate q^2 to get b-value
        self.message("{:.3f}  {:.2f}  {:.2f}  {:.2f}".format(Imav, Imav*self.GCal*1000, q180, bvalue))
        self.measuredBValues[i]=bvalue
      if self.DACmax!= 0:
        if self.gradOrientation=='X':      #Set gradient orientation and max gradient current when DAC is set to 100
            self.Axmax=Imav*100/self.DACmax
            self.setGradOrientation('X')
            self.message("Reset Gx Imax={:.3f}A".format(self.Amax))
        if self.gradOrientation=='Y':
            self.Aymax=Imav*100/self.DACmax
            self.setGradOrientation('Y')
            self.message("Reset Gy Imax={:.3f}A".format(self.Amax))
        if self.gradOrientation=='Z':
            self.Azmax=Imav*100/self.DACmax
            self.setGradOrientation('Z')
            self.message("Reset Gz Imax={:.3f}A".format(self.Amax))                     
      self.dataImg.clear()
      self.dataReal.clear()
      self.dataReal.setTitle('Gradient Waveforms' , color='b', size='18pt')
      self.dataImg.setTitle('q/2&pi; (1/mm)' , color='b', size='18pt')
      self.dataReal.setLabel('bottom','Time(ms)' , **self.labelStyle)
      self.dataImg.setLabel('bottom', 'Time(ms)', **self.labelStyle)
      self.dataReal.setLabel('left','Gradient(T/m)' , **self.labelStyle)
      self.dataImg.setLabel('left', 'q/2&pi;(1/mm)', **self.labelStyle)
      self.dataReal.plotItem.getAxis('left').setPen(self.wpen)
      self.dataReal.plotItem.getAxis('bottom').setPen(self.wpen)
      self.dataImg.plotItem.getAxis('left').setPen(self.wpen)
      self.dataImg.plotItem.getAxis('bottom').setPen(self.wpen)            
      for i, ctrace in enumerate(self.gradCurrentTraces):
        p=pg.mkPen(pg.intColor(i+self.penstep), width=2)
        self.dataReal.plot(ctrace[:,0], ctrace[:,1], pen=p, width=2, name='current trace')   #plot gradient traces
      for i, ctrace in enumerate(self.gradQTraces):
        p=pg.mkPen(pg.intColor(i+self.penstep), width=2)
        self.dataImg.plot(ctrace[0,:],ctrace[1,:], pen=p, width=2, name='current trace')   #plot wavevector q/2pi
      self.addDataPlottoReport()        
      self.addCrossHairs()


  def changeGradOrientation(self):
      direction=self.ui.cbGradOrientation.currentText()
      self.setGradOrientation(direction)
      
  def setGradOrientation(self, direction):
        if direction=='X':
            index = self.ui.cbGradOrientation.findText('X', Qt.MatchFixedString)     #set data type in gui combo box
            self.gradOrientation='X'
            self.Amax=self.Axmax
            self.GCal=self.GxCal
            self.ui.leGradCal.setText('{:10.3f}'.format(1000*self.GCal))     #grad calibrations are in T/m/a, display is in mT/m/A
            self.ui.leImax.setText("{:.3f}".format(self.Axmax))
        if direction=='Y':
            index = self.ui.cbGradOrientation.findText('Y', Qt.MatchFixedString)     #set data type in gui combo box
            self.gradOrientation='Y'
            self.Amax=self.Aymax
            self.GCal=self.GyCal
            self.ui.leGradCal.setText('{:10.3f}'.format(1000*self.GCal))
            self.ui.leImax.setText("{:.3f}".format(self.Aymax))    
        if direction=='Z':
            index = self.ui.cbGradOrientation.findText('Z', Qt.MatchFixedString)     #set data type in gui combo box
            self.gradOrientation='Z'
            self.Amax=self.Azmax
            self.GCal=self.GzCal
            self.ui.leGradCal.setText('{:10.3f}'.format(1000*self.GCal))
            self.ui.leImax.setText("{:.3f}".format(self.Azmax))
        if index >= 0:
            self.ui.cbGradOrientation.setCurrentIndex(index)
            
  def calculateImageWidths(self):
    '''Calculates the frequency width df of 1d images of object of width d mm, 
    to use in gradient G calibrations df=gamma G d/2pi, G=gcal I'''
    #QMessageBox.about(self, "1D Image Width and Gradient Calibration", "1. Calculates frequency width, df, of 1d images of object of width d,<br> 2. Calibrates gradients using: df=gamma G d/2pi, G=gcal I (G=gradient, I=current, gcal=grad cal factor)")
    self.message('<b>Calculate Image Widths: </b> ' + self.gradOrientation + '-direction')
    self.openDataFile(message='Open 1D Image for Gradient Calibration ' +self.gradOrientation + '-direction')
    self.gradientDirectory=self.direc     #Save gradient directory name
    self.subtractBaselines()
    self.FFTData()
    try:
      dat=np.absolute(self.data[:,:,0,0])     #fit magnitude of FT data
    except:
        QMessageBox.about(self,'1D image data not found',"Load data")
        return None  
    self.gradCal=np.zeros(self.data.shape[1])   #array of gradient values for each gradient current
    self.gradCalerror=np.zeros(self.data.shape[1])   #array of gradient error values taken form 1d image fit sandard error for each gradient current
    f=self.freq
    res=np.zeros((self.data.shape[0],self.data.shape[1]))    #array to store residuals from fit
#     nf=np.int(dat.shape[0]/4)        #use threshhold to determine width
#     th=np.max(dat[0:nf])*5
#     self.message('threshold={:.3e}'.format(th))
#     self.inf2.setValue(th)
    dx, OK= QInputDialog.getDouble(self,"Input cell dimension", "cell diameter: dxGCal(mm)", value=self.dxgcal*1000, decimals=3)
    if OK:
        self.dxgcal=dx/1000
    self.message('<b>Cell diameter: </b>' +self.gradOrientation+ '-direction d(mm)={:.3f}'.format(self.dxgcal*1000))
    self.message('Igrad(A), Smax,  f0(Hz),  df(Hz), bg, G(mT/m),  Gerr(mT/m), ferr (Hz), DAC')
    if self.gradOrientation=='Z':
        ax=0
    else:
        ax=1
    for i in range(self.data.shape[1]):
        params=image1dfit.initialize (f=f, data=dat[:,i],ax=ax)
        pdicti=params[0] #parameter dictionary
        plist=params[1] #parameter list
        fitoutput = lmfit.minimize(image1dfit.ObjFnc,pdicti,args=(f,dat[:,i]))
        pdict=fitoutput.params
        A=pdict['A'].value
        f0=pdict['f0'].value
        gf=pdict['gf'].value
        gfStdEr=pdict['gf'].stderr
        bl=pdict['bl'].value
        res[:,i]=fitoutput.residual
#       aboveThresh=np.argwhere(dat[:,i]>th)
#       minf=self.freq[np.min(aboveThresh)]
#       maxf=self.freq[np.max(aboveThresh)]
#       df=np.absolute(maxf-minf)
        self.gradCal[i]=np.sign(self.measuredI[i])*2*np.pi*2*gf/self.Gamma/self.dxgcal  #calculate gradient from df
        self.gradCalerror[i]=2*np.pi*2*gfStdEr/self.Gamma/self.dxgcal  #calculate gradient fit error
        DAC=50*self.GrAmp*self.gradAmpArray[i] #cacluate maximum DAC output for gradient pulse
        self.message('{:10.3f} {:6.1e} {:10.1f} {:10.1f} {:10.1f}  {:10.2f}  {:10.2f}  {:10.2f}  {:10.0f}'.format(self.measuredI[i],A,f0,2*gf,bl, self.gradCal[i]*1000, self.gradCalerror[i]*1000, gfStdEr, DAC))
    self.GCal, self.GcalOffset=np.polyfit(self.measuredI,self.gradCal,1)      #fit line to gradient vs current to get Gcal
    self.Amax=abs(self.measuredI[-1]-self.measuredI[0])/2/(self.GrAmp/100)        #Amax is the (maximum positive current - max negative current)/2 when GrAmp=100,  
    if self.gradOrientation=='X':   #update gradient calibrations
        self.Axmax=self.Amax
        self.GxCal=self.GCal
    if self.gradOrientation=='Y':
        self.Aymax=self.Amax
        self.GyCal=self.GCal
    if self.gradOrientation=='Z':
        self.Azmax=self.Amax
        self.GzCal=self.GCal
    self.message('<b>Set:</b> gradient direction=' + self.gradOrientation + ', Gcal(mT/m/A)=' + '{:10.3f}'.format(1000*self.GCal) +', Imax(A)=' + '{:10.3f}'.format(self.Amax))
    self.ui.leGradCal.setText('{:10.3f}'.format(1000*self.GCal))
    self.ui.leImax.setText('{:10.4f}'.format(self.Amax))
 
    self.image1DPlot=self.plotWindow(self)      #open a plot window for 1dimages and fits
    self.image1DPlot.win.setWindowTitle('1D Image ' + self.fileName)
    self.image1DPlot.dplot.setLogMode(False, True)
    for i in range(self.nSpectra):
              p=pg.mkPen(pg.intColor(i+self.penstep), width=2)
              self.image1DPlot.dplot.plot(f,dat[:,i], pen=p)
              self.image1DPlot.addData(f, dat[:,i], description='1dFit Images')
              self.image1DPlot.addResData(f, res[:,i], description='1dFit Residuals')
    self.image1DPlot.dplot.setLabel('bottom', "Frequency (Hz)") #**self.labelStyle)
    self.image1DPlot.dplot.setLabel('left', 'Signal')#, **self.labelStyle)
    self.image1DPlot.dplot.setTitle('1d Image in '+ self.gradOrientation,size='14pt', color=(100,200,200))
    self.image1DPlot.show() 
    self.image1DPlot.exportPlot()       #make a copy in the report
    
    
    self.gcalPlot=self.plotWindow(self)     #make a plot of gradient strenght vs current along with fit
    self.gcalPlot.win.setWindowTitle('Gradient Calibration ' + self.fileName)
    self.gcalPlot.dplot.setLogMode(False, False)
    self.gcalPlot.dplot.plot(self.measuredI,self.gradCal, pen=None, symbol='o',symbolSize=14,)
    self.gcalPlot.dplot.plot(self.measuredI,self.measuredI*self.GCal+self.GcalOffset)
    self.gcalPlot.dplot.setLabel('bottom', "Current (A)") #**self.labelStyle)
    self.gcalPlot.dplot.setLabel('left', "Gradient (T/m)") #, **self.labelStyle)
    self.gcalPlot.plotTitle='{} Gradient Calibration: Gcal(mT/m/A)={:10.3f}'.format(self.gradOrientation,1000*self.GCal )
    self.gcalPlot.setPlotTitle()
    self.gcalPlot.show() 
    self.gcalPlot.exportPlot()  #make a copy in the report

  def inputPGSEKurtosis(self):
      self.inputPGSEDiffusion(kurtosis=True)
  
  def inputPGSEDiffusion(self,kurtosis=False, biExp=False):
    ''' opens NMR diffusion data file, processes it using standard method'''
    if self.ui.cbDiffusionModel.currentText()=='Bi-exponential':
      biExp=True
    elif self.ui.cbDiffusionModel.currentText()=='Kurtosis':
      kurtosis=True
    self.openDataFile(message='Open Diffusion PGSE File, Gradient direction= '+self.gradOrientation)
    self.subtractBaselines()
    self.FFTData()
    self.setPhase(firstSIndex=True)     #find maximum and phase to set imag to zero at that point
    self.findIntWidth() #needed to set the range for maximized integral of real part of spectra
    self.finePhase()
    self.hideCrossHairs()
    self.addDataPlottoReport()
    self.integrateData()
    self.fitDiffusion(Kurtosis=kurtosis,biExp=biExp)
    self.pw.exportPlot()    #puts a copy of the plot into the report
    
  def inputEddyCurrentCorrection(self):
    ''' opens NMR diffusion data file, processes it using standard method and stors result as an eddy current correction'''
    self.openDataFile(message='Open Eddy Current Correction PGSE File, Gradient direction= '+self.gradOrientation)
    self.subtractBaselines()
    self.FFTData()
    self.setPhase()
    self.findIntWidth()
    self.finePhase()
    self.addDataPlottoReport()
    self.integrateData()
    self.pw.fitDiffusion()
    self.pw.exportPlot()
    self.updateECC()

  def updateECC(self):
    direction=self.ui.cbGradOrientation.currentText()    
    if direction=='X':
      self.eddyGxCorrection=self.ADCreported      #cupdate current ECC with value measured from ECC pulse sequence
    if direction=='Y':
      self.eddyGyCorrection=self.ADCreported      #cupdate current ECC with value measured from ECC pulse sequence
    if direction=='Z':
      self.eddyGzCorrection=self.ADCreported      #cupdate current ECC with value measured from ECC pulse sequence
    self.ui.leECCX.setText("{:.4f}".format(self.eddyGxCorrection*1000))
    self.ui.leECCY.setText("{:.4f}".format(self.eddyGyCorrection*1000))          


  def findIntWidth(self):
    '''Find the integration width about the bigget peak in real part of the spectra using the selected number of peak widths self.npws or a fixed bandwidth '''
    self.dataType=self.ui.cbPulseSequence.currentText()
    fl=0
    if self.dataType=='IR':     #pick the last curve longest TI to find peak width
      fl=-1
    if self.dataType=='CPMG':   #pick first curve with the shorted aquisition time to find peak width
      fl=0
    if self.dataType=='Diffusion':   #pick first curve with the smallest bValue time to find peak width
      fl=0
    if self.dataType=='nutation'or self.dataType=='RFattn':       #Find the spectra with the largest value, use this spectra to find integration range
      smax=np.amax(self.data[:,:,self.repeatIndex,0].real, axis=0)
      fl=np.argmax(smax)
    self.npws=int(self.ui.lenPW50.text())   #number of fullwidth at 50% max used for integration range
    y=self.data[:,fl,self.repeatIndex,0].real
    maxVal=np.amax(y)   #find the maximum value, will select region around this maximum value to integrate
    maxInd=np.argmax(y)
    maxVal50 = 0.5*maxVal
    biggerCondition = [a > maxVal50 for a in y] #returns boolean array which is true for values above 50% of max value
    width=np.sum(biggerCondition)   #returns nuber of points above 50 max  ***Needs to be fixed if there are several large peaks***
    df=1/self.tntfile.dwell[0]/self.data.shape[0]   #frequency spacing in Hz
    self.FirstPeakfwhm=width*df
    self.T2star=1/(np.pi*self.FirstPeakfwhm)
    self.message('<b>Integration:</b> from spectra={}, peak found at {},f0={:.4f}Hz, approx FWHM={:.3f}Hz, Smax={:.3e}, T2*(ms)={:.5f}'.format(fl,maxInd,self.freq[maxInd],self.FirstPeakfwhm,self.FirstPeakfwhm/self.Gammaf,self.T2star*1000))
    if self.ui.cbIntegrationRange.currentText()=='Integrate n linewidths':
        wind=self.npws*width
    if self.ui.cbIntegrationRange.currentText()=='Integrate fixed bandwidth':
        self.intbw=float(self.ui.leIntBW.text())
        wind=int(self.intbw/df/2)    
    self.message(self.ui.cbIntegrationRange.currentText()+', <b>Integration width (Hz)</b>={:.5f}'.format(2*wind*df))
    minI=maxInd-wind
    if minI<0:
        minI=0
    maxI=maxInd+wind
    if maxI>self.data.shape[0]-1:
        maxI=self.data.shape[0]-1
    self.ui.leiStart.setText(str(minI))
    self.ui.lejStop.setText(str(maxI))
    self.ui.leFstart.setText('{:6.3f}'.format(self.freq[minI]))
    self.ui.leFstop.setText('{:6.3f}'.format(self.freq[maxI]))



  def setPhase(self, firstSIndex=False):
    '''Set phase for all spectra using simplest algorythym find maximum signal phase to make real part maximum and imaginary zero'''
    ls=int(self.data.shape[0]/20)
    if firstSIndex:     #reduce the peak search area to a region around the first peak maxima,
        smax= np.argmax(np.abs(self.data[:,0,0,0]))
        il=smax-ls
        if il<0:
            il=0
        ih=smax+ls
        if ih>self.data.shape[0]-1:
            ih=self.data.shape[0]-1
    else:
        il=0
        ih=self.data.shape[0]-1
    for j in range(self.nRepeats):
        for i in range(self.nSpectra):
            marg=il+np.argmax(np.abs(self.data[il:ih,i,j,0]))
            dphase=np.angle(self.data[marg,i,j,0])
            self.data[:,i,j,0]= self.data[:,i,j,0]*np.exp(-1j*dphase)
            self.PhaseArray[i,j]=180*dphase/np.pi
            self.FreqArray[i,j]=self.freq[marg]
    #self.Phase=None
    self.plotData()  #plot data
    self.message('<b>Phase angles=</b>'+ np.array2string(self.PhaseArray[:,0], precision=2, separator='  ',suppress_small=True).strip("\n\t "))
    self.message('<b>Peak frequency=</b>'+ np.array2string(self.FreqArray[:,0], precision=3, separator='  ',suppress_small=True).strip("\n\t "))


  def finePhase(self):
    '''fine phase adjust for all spectra, maximizes integral of real part of spectrum'''
    self.iStart=int(self.ui.leiStart.text())
    self.jStop=int(self.ui.lejStop.text())
    for j in range(self.nRepeats):
        for i in range(self.nSpectra):
            anglemax=self.finePhaseIJ(i,j)
            self.data[:,i,j,0]= self.data[:,i,j,0]*np.exp(-1j*anglemax*np.pi/180)
            self.PhaseArray[i,j]+=anglemax
    self.plotData()  #plot data
    for j in range(self.nRepeats):
        self.message('<b>Fine phase adjust {}: Phase angles=</b>'.format(j)+ np.array2string(self.PhaseArray[:,j], precision=2, separator=' ',suppress_small=True).strip("\n\t "))

  def finePhaseIJ(self, i=0,j=0, maxparam='ssum', angleRange=45):
    '''finds and returns the phase angle that maximizes the integral of the real part of spectra IJ'''
    smax=0.0
    anglemax=0.0
    angle=np.arange(-angleRange,angleRange)
    ssum=np.zeros(len(angle))
    nsum=np.zeros(len(angle))
    psum=np.zeros(len(angle))
    for ind, a in enumerate(angle):
        phased=self.data[self.iStart:self.jStop,i,j,0]*np.exp(-1j*a*np.pi/180)      #phase shift data, see if the real part is more positive
        ssum[ind]=abs(np.sum(phased.real))
        nsum[ind]=((phased.real<0)*phased.real).sum()
        psum[ind]=((phased.real>0)*phased.real).sum()
    if maxparam=='ssum':    #return angle which gives max of integral
        anglemax=angle[np.argmax(ssum)]            
    if maxparam=='nsum':    #return angle which is  min of negative component, make a postive real spectra nonnegative
        anglemax=angle[np.argmax(nsum)]
    if maxparam=='psum':    #return angle which is  min of positive compnent, make a neagtive real spectra nonpositive eg first spectra in T1IR
        anglemax=angle[np.argmin(psum)]
    return anglemax

  def nutation(self):
      '''opens nutation file, baseline subtract, FFT, phasing, int width, and plot signal vs RF pulse time'''
      self.message('******Nutation Processing******', bold=True, color='blue')
      s=self.openDataFile(message='Open nutation file')
      if s=='cancel':
          return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse phase adjust
        i=5 # Phase on first spectra
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase) #coarse phase so ith spectra real part is maxima and imaginary part close to zero
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='psum') #fine phase adjust, make real part most negative 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase 
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.integrateData()
      self.fitNutation()
      self.pw.exportPlot() 

  def noiseAnalysis(self):
      '''opens FID file, expects a single FID in dimension 1, baseline subtract, Opens new noise plot with FID phase data, '''
      self.message('******Noise Processing******', bold=True, color='blue')
      s=self.openDataFile(message='Open FID file for noise analysis')
      if s=='cancel':
          return
      self.subtractBaselines()
      self.noiseData=copy.copy(self.data[:,0,0,0])
      self.ftNoiseData=fft.fftshift(fft.fft(self.noiseData))       #find freqencey of largest pek, shift that frequency to 0
      smax= np.argmax(np.abs(self.ftNoiseData))
      f=self.freq[smax]     #Frequency at peak maximum
      self.message('<b>Center Noise Spectra</b> by shifting by {:4.2f}Hz'.format(f))
      self.noiseData= self.noiseData*np.exp(1j*2*np.pi*f*self.tfid) #multiply noise FID to get to zero field
      self.data[:,0,0,0]=self.data[:,0,0,0]*np.exp(1j*2*np.pi*f*self.tfid)
      self.plotData(showCrossHairs=False)  #plot data
      self.noiseWindow=self.plotNoise(self)
      self.noiseWindow.win.setWindowTitle('Noise Data' + self.fileName)
      self.noiseWindow.show()
      #self.noiseWindow.setLogMode(logx=False,logy=False)
      self.noiseWindow.dplot.setLabel('bottom', "Time") 
      self.noiseWindow.magData=np.absolute(self.noiseData)
      self.noiseWindow.phaseData=np.angle(self.noiseData) 
      self.ftNoiseData=fft.fftshift(fft.fft(np.unwrap(self.noiseWindow.phaseData))) 
      self.noiseWindow.tfid=self.tfid
      self.noiseWindow.maxIndex=len(self.tfid)
      self.noiseWindow.freq=self.freq
      self.noiseWindow.dwellTime=self.dwellTime
      self.noiseWindow.setMaxTime(tMax=2)       #sets upper time of interest to 2s and plots phase
      self.RIorAP='AP'  #set display to amplitude/phase
      self.ui.sbViewCurveN.setValue(0)
      self.dataReal.setLogMode(False,True)
      self.FFTData()    #FFTs data to get a spectrum and plots spectrum


            
  def  ECCGradientRingdown(self):
      '''opens  ECCGradientRingdown file, baseline subtract, FFT, phasing, int width, and plots resoant frequency and inewidth vs ringdown time'''
      self.message('******ECCGradientRingdown Processing******', bold=True, color='blue')
      s=self.openDataFile(message='Open ECCGradientRingdown file')
      if s=='cancel':
          return
      self.subtractBaselines()
      self.FFTData()
      self.setPhase()
      self.plotData(showCrossHairs=False)  #plot data
      self.findIntWidth()
      self.integrateData()
      #self.fitNutation()
      #self.pw.exportPlot() 
        
  def T1IR(self):
      '''opens T1IR file, baseline subtract, FFT, phasing, int width, and plot signal vs inversion time'''
      self.message('******T1 Inversion Recovery Processing******', bold=True, color='blue')
      s=self.openDataFile(message='Open T1IR file')
      if s=='cancel':
          return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse phase adjust
        i=0 # Phase on first spectra
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= -self.data[:,:,j,0]*np.exp(-1j*dphase) #coarse phase so first spectra real part is maxima and imaginary part close to zero
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='psum') #fine phase adjust, make real part most negative 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase 
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.integrateData()
      self.pw.T1Fit()
      self.pw.exportPlot()
      
  def T1IRSpectra(self):
      '''opens T1IR file, baseline subtract, FFT, phasing, int width, and plot signal vs inversion time'''
      self.message('******T1 Inversion Recovery Processing******', bold=True, color='blue')
      s=self.openDataFile(message='Open T1IR file')
      if s=='cancel':
          return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse phase adjust
        i=0 # Phase on first spectra
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= -self.data[:,:,j,0]*np.exp(-1j*dphase) #coarse phase so first spectra real part is maxima and imaginary part close to zero
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='psum') #fine phase adjust, make real part most negative 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase 
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      m=4
      self.reddata=self.data[:,:,:,0].real.mean(axis=2)#define reduced data as a 2 dimensional array averaged over measurements
      self.reddata/=np.amax(self.reddata)     #normalize reduced data
      pc=[]     #pointcloud for T1, T2 spectra, contains frequency, component amplitude, T1
      for i in range (2,self.nPoints-2):
          data=(self.reddata[i-2,:]+self.reddata[i-1,:]+self.reddata[i,:]+self.reddata[i+1,:]+self.reddata[i+2,:])/5
          if np.amax(data)> self.spectraMinimum:
              params=T1IR.initializeT1IR (TI=self.TI, data=data)
              pdicti=params[0] #parameter dictionary
              plist=params[1] #parameter list
              fitoutput = lmfit.minimize(T1IR.T1IR,pdicti,args=(self.TI ,data))
              pdict=fitoutput.params
              if pdict['T1'].value<10 and pdict['A'].value>0.01:        #Accept points in which the relaxation time is in range 0 to 10 s and the normalized amplitude is greater than 1%
                  p=np.array([self.freq[i],pdict['A'].value,1000*pdict['T1'].value])
                  pc.append(p)
      pc=np.array(pc)
      cmax=pc[:,1].sum()
      pc[:1]/=cmax  #normalize exponetial amplitudes so integral = 1
      npoints=pc.shape[0]
      colorpoints=255*np.column_stack((pc[:,1],pc[:,1],pc[:,1]))#,np.full(npoints,0.5)))     #4 column color array
      color=np.column_stack((pc[:,1],pc[:,1],pc[:,1],np.full(npoints,0.5)))     #4 column color array
      self.plotPointCloud(pc, color=color)
      self.pw=self.plotWindow(self)
#      cm = pg.colormap.get('CET-L3') # prepare a linear color map
      self.pw.win.setWindowTitle('T1 Spectra :' + self.fileName)
      self.pw.dplot.setLabel('bottom', "Frequency(Hz)")
      self.pw.dplot.setLabel('left', "T1(ms)")
      self.pw.show()
      self.pw.dplot.plot(pc[:,0],pc[:,2], pen=None, symbol='o',symbolSize=10, symbolBrush =[pg.mkBrush(v) for v in pc[:,1]])  #symbolBrush=(0,0,225)
      self.pw.dplot.plot(pc[:,0],pc[:,1], pen=None, symbol='o',symbolSize=0, symbolBrush =(50,0,225))  #symbolBrush=(0,0,225)
      self.pw.dplot.plot(self.freq,-1000*self.data[:,0,0,0].real/np.amax(-self.data[:,0,0,0].real), pen=self.bpen) #plot spectra for reference

  def T1IRbiBW(self):
      '''Fits T1 data with narrow, then wide bandwidth to give major peak T1 and effective MRI T1
      opens T1IR file, baseline subtract, FFT, phasing, int widt, and plot signal vs inversion time'''
      self.message('******T1 Inversion Recovery Processing******', bold=True, color='blue')
      s=self.openDataFile(message='Open T1IR file')
      if s=='cancel':
          return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse phase adjust
        i=0 # Phase on first spectra
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= -self.data[:,:,j,0]*np.exp(-1j*dphase) #coarse phase so first spectra real part is maxima and imaginary part close to zero
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='psum') #fine phase adjust, make real part most negative 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase 
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.ui.cbIntegrationRange.setCurrentIndex(self.ui.cbIntegrationRange.findText('Integrate fixed bandwidth'))      #set to integrate over fixed bandwith, first narrow, then wide
      self.intbw=self.intbwNarrow
      self.findIntWidth()
      self.integrateData()
      self.pw.T1Fit()
      self.pw.exportPlot()
      T1Narrow=(self.T1reported*1000, self.T1StdErrReported*1000, self.deltaReported,self.intbw) 
      self.intbw=self.intbwWide
      self.findIntWidth()
      self.integrateData()
      self.pw.T1Fit()
      self.pw.exportPlot()
      T1Wide=(self.T1reported*1000, self.T1StdErrReported*1000, self.deltaReported,self.intbw)
      self.message('<b>Narrow BW:</b> T1(ms), T1 stder(ms),delta, bw(Hz)={:6.1f} {:6.1f} {:6.3f} {:6.1f}'.format(T1Narrow[0],T1Narrow[1],T1Narrow[2],T1Narrow[3]))
      self.message('<b>Wide BW:</b> T1(ms), delta, bw(Hz)={:6.1f} {:6.1f} {:6.3f} {:6.1f}'.format(T1Wide[0],T1Wide[1],T1Wide[2],T1Wide[3]))
                        
  def T2CPMG(self):
      '''opens T2CPMG file, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******T2 CPMG Processing******', bold=True, color='blue')
      c=self.openDataFile(message='Open T2CPMG file')
      if c=='cancel':
        return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use first spectra  
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.integrateData()
      self.pw.parWin.fitCPMGT2()
      self.pw.exportPlot()

  def T2CPMG_tau_array(self):
      '''opens T2CPMG tau_array file, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******T2 CPMG Processing******', bold=True, color='blue')
      c=self.openDataFile(message='Open T2CPMG file')
      if c=='cancel':
        return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use first spectra  
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax
        s+='tau{}(ms)={:4.1f}, tau phase(deg)= {:4.1f};  '.format(j,1000*self.cpmg_tau_array[j],self.PhaseArray[i,j]) 
      self.message(s)
      self.message('')
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.integrateData()
      self.pw.parWin.fitCPMGT2_tau_array()
      self.T2Window=self.generalPlot(self)
      self.T2Window.win.setWindowTitle('CPMG T2 vs refocus freqeuncy' + self.fileName)
      self.T2Window.show()
      #self.noiseWindow.setLogMode(logx=False,logy=False)
      self.T2Window.dplot.plot(0.5/(self.cpmg_tau_array+self.tauCPMG), self.T2fit*1000.0, symbol='o')
#      self.pw.exportPlot() 
           
  def T2CPMGbiExp(self):
      '''opens T2CPMG file for biExponential fit, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******T2 CPMG Processing******', bold=True, color='blue')
      c=self.openDataFile(message='Open T2CPMG file')
      if c=='cancel':
        return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive 
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.integrateData()
      self.fitCPMGT2biEx()
      self.pw.exportPlot()
      
  def T2CPMGbiExpNWBW(self):
      '''opens T2CPMG file fits to biExponetial using both narrow and wide bandwidths, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******T2 CPMG biExponential Processing******', bold=True, color='blue')
      c=self.openDataFile(message='Open T2CPMG file')
      if c=='cancel':
        return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive 
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.ui.cbIntegrationRange.setCurrentIndex(self.ui.cbIntegrationRange.findText('Integrate fixed bandwidth'))      #set to integrate over fixed bandwith, first narrow, then wide
      self.intbw=self.intbwNarrow
      self.ui.leIntBW.setText('{:.1f}'.format(self.intbw))
      self.findIntWidth()
      self.integrateData()
      self.fitCPMGT2biEx()
      bwtext=pg.TextItem(text = 'NarrowBW(Hz)={:.1f}'.format(self.intbw), color=(150, 200, 100))
      self.pw.dplot.addItem(bwtext)
      bwtext.setPos(.3, 0)
      self.pw.exportPlot()
      T2Narrow=(self.T2areported*1000, self.T2aStdErrReported*1000, self.CaReported*100,self.T2breported*1000,self.T2bStdErrReported*1000, self.CbReported*100,self.intbw) 
      self.intbw=self.intbwWide
      self.ui.leIntBW.setText('{:.1f}'.format(self.intbw))
      self.findIntWidth()
      self.integrateData()
      self.fitCPMGT2biEx()
      bwtext=pg.TextItem(text = 'WideBW(Hz)={:.1f}'.format(self.intbw), color=(150, 200, 100))
      self.pw.dplot.addItem(bwtext)
      bwtext.setPos(.3, 0)
      self.pw.exportPlot()
      T2Wide=(self.T2areported*1000,self.T2aStdErrReported*1000, self.CaReported*100,self.T2breported*1000,self.T2bStdErrReported*1000, self.CbReported*100,self.intbw) 
      self.message('<b>Narrow BW:</b> T2a(ms)  T2aStdErr(ms)  a-fraction(%)  T2b(ms)  T2abStdErr(ms)  b-fraction(%)  bw(Hz)={:6.1f} {:6.1f} {:6.1f} {:6.1f} {:6.1f} {:6.1f} {:6.1f}'.format(T2Narrow[0],T2Narrow[1],T2Narrow[2],T2Narrow[3],T2Narrow[4],T2Narrow[5],T2Narrow[6]))
      self.message('<b>Wide BW:</b> T2a(ms)  T2aStdErr(ms)  a-fraction(%)  T2b(ms) T2bStdErr(ms)  b-fraction(%)  bw(Hz)={:6.1f} {:6.1f} {:6.1f} {:6.1f} {:6.1f} {:6.1f} {:6.1f}'.format(T2Wide[0],T2Wide[1],T2Wide[2],T2Wide[3],T2Wide[4],T2Wide[5],T2Wide[6]))
                        
  def T2Spectra(self, spectralAveWidth=5):
      '''opens T2CPMG file, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******T2 Spectra Processing******', bold=True, color='blue')
      self.openDataFile(message='Open T2CPMG file')
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive       
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax 
        self.message('Repeat {}, total phase adjust (deg)={:4.1f}'.format(j,self.PhaseArray[i,j]))
      self.plotData()  #plot data
      m=4
      #self.reddata=self.data[:,:,:,0].real.mean(axis=2)#define reduced data as a 2 dimensional array averaged over measurements
      self.reddata=self.data[:,:,0,0].real
      self.reddata/=np.amax(self.reddata)     #normalize reduced data
      self.tauCPMG=float(self.ui.leTauCPMG.text())*1E-3     #read in tauCPMG in ms convert to aquition time in s
      self.cpmgT180=float(self.ui.leT180.text())*1E-6
      self.taCPMG=(2.0*(self.tauCPMG+self.cpmgDelay) +self.cpmgT180)*self.cpmgLoopIt
      pc=[]     #pointcloud for T1, T2 spectra, contains frequency, Si compnent amplitude, T2
      for i in range (spectralAveWidth,self.nPoints-spectralAveWidth):
          #data=(self.reddata[i-2,:]+self.reddata[i-1,:]+self.reddata[i,:]+self.reddata[i+1,:]+self.reddata[i+2,:])/5
          data=np.sum(self.reddata[i-spectralAveWidth:i+spectralAveWidth+1,:], axis=0)/(2*spectralAveWidth+1)
          if np.amax(data)> self.spectraMinimum:
              params=T2CPMGbiExp.initialize (TE=self.taCPMG, data=data)
              pdicti=params[0] #parameter dictionary
              plist=params[1] #parameter list
              fitoutput = lmfit.minimize(T2CPMGbiExp.ObjFnc,pdicti,args=(self.taCPMG,data))
              pdict=fitoutput.params
              if pdict['T2a'].value<10 and pdict['Si'].value>0.01:
                  p=np.array([self.freq[i],pdict['Si'].value,1000*pdict['T2a'].value])
                  pc.append(p)
              if pdict['T2b'].value<10 and pdict['Si'].value>0.01:
                  p=np.array([self.freq[i],pdict['Sib'].value,1000*pdict['T2b'].value])
                  pc.append(p)
      pc=np.array(pc)
      cmax=pc[:,1].sum()
      pc[:1]/=cmax  #normalize exponential amplitudes so integral = 1
      npoints=pc.shape[0]
      colorpoints=255*np.column_stack((pc[:,1],pc[:,1],pc[:,1]))#,np.full(npoints,0.5)))     #4 column color array
      color=np.column_stack((pc[:,1],pc[:,1],pc[:,1],np.full(npoints,0.5)))     #4 column color array
      #self.plotPointCloud(pc, color=color)
      self.pw=self.plotWindow(self)
#      cm = pg.colormap.get('CET-L3') # prepare a linear color map
      self.pw.win.setWindowTitle('T2 Spectra :' + self.fileName)
      self.pw.dplot.setLabel('bottom', "Frequency(Hz)")
      self.pw.dplot.setLabel('left', "T2(ms)")
      self.pw.show()
      self.pw.dplot.plot(pc[:,0],pc[:,2], pen=None, symbol='o',symbolSize=10, symbolBrush =(50,0,225))  #symbolBrush=(0,0,225)
      self.pw.dplot.plot(pc[:,0],pc[:,1], pen=None, symbol='o',symbolSize=0, symbolBrush =(50,0,225))  #symbolBrush=(0,0,225)
      self.pw.dplot.plot(self.freq,200*self.data[:,0,0,0].real/np.amax(self.data[:,0,0,0].real), pen=self.bpen) #plot spectra for reference

      
  def T2SE(self):
      '''opens T2SEfile, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******T2SE Processing******', bold=True, color='blue')
      self.openDataFile(message='Open T2SE file')
      self.subtractBaselines()
      self.FFTData()
      self.setPhase()
      # for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
      #   i=0
      #   marg=np.argmax(np.abs(self.data[:,i,j,0]))
      #   dphase=np.angle(self.data[marg,i,j,0])
      #   self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
      #   self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive       
      # for j in range(self.nRepeats):    #fine phase adjust   
      #   self.repeatIndex=j
      #   anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
      #   self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
      #   self.PhaseArray[:,j]+=anglemax 
      #   self.message('Repeat {}, total phase adjust (deg)={:4.1f}'.format(j,self.PhaseArray[i,j]))
      self.plotData(showCrossHairs=False)  #plot data
      self.integrateData()
      self.pw.parWin.fitCPMGT2()
      
  def T1rho(self):
      '''opens T1rho file, baseline subtract, FFT, phasing, int widt, and plot signal vs spinlock time time'''
      self.message('******T1rho Processing******', bold=True, color='blue')
      self.openDataFile(message='Open T1rho file')
      self.subtractBaselines()
      self.FFTData()
      self.setPhase()
      self.findIntWidth()
      self.finePhase()
      self.plotData(showCrossHairs=False)  #plot data
      self.integrateData()
      self.fitT1rho()

  def setDiffusionAnalysisType(self):
      if self.ui.chGradCal.isChecked():    #Flag to determine if gradient calibratison and ECC should be done, otherwise use default grad calibrations and ECC=0
        self.performGradientCal = True
      else:
        self.performGradientCal = False
      if self.ui.chECCCal.isChecked():    #Flag to determine if gradient calibratison and ECC should be done, otherwise use default grad calibrations and ECC=0
        self.performECC = True
      else:
        self.performECC = False
      
  def gradCalXY(self):
    '''Calibrate XY gradients'''
    self.message('******Gradient Calibration******', bold=True, color='blue')
    self.setDiffusionAnalysisType()
    try:
        for grador in ['X','Y']:
            self.setGradOrientation(grador)
            self.inputGradientCurrentTraces()
            self.calculateImageWidths()
        self.message('ImaxGx, ImaxGy(A)={:.3f} {:.3f}, Gradient Cal: Gxcal, Gycal(mT/m/A)={:.3f} {:.3f}'.format(self.Axmax,self.Aymax,1000*self.GxCal, 1000*self.GyCal), color='green', bold=True)
        txtgradcal='ImaxGx(A)={:.3f}, ImaxGy(A)={:.3f}, Gxcal(mT/m/A)={:.3f}, Gycal(mT/m/A)={:.3f} \n'.format(self.Axmax,self.Aymax,1000*self.GxCal, 1000*self.GyCal)
        fname=self.gradientDirectory + '\\GradientCal.txt'
        f=open(fname, 'w')
        f.write(txtgradcal)
        f.close()
        self.message('Gradient calibrations written to ' + fname +'\n')
    except:
      self.message('Gradient Calibration Error')
      raise

      
  def diffusionAnalysis(self, biExp=False):
    '''Diffusion analysis: inputs gradient currents and 1d images to calibrate gradients, then eddy current correction data, and finaly diffusion data'''
    self.message('******Diffusion Processing******', bold=True, color='blue')
    self.doNotFitDBelowSNR=True
    self.setDiffusionAnalysisType()
    try:
        for grador in ['X','Y']:
            self.setGradOrientation(grador)
            if self.performGradientCal == True:
                self.inputGradientCurrentTraces()
                self.calculateImageWidths()
                
            if self.performECC == True:
                self.inputEddyCurrentCorrection()
            self.inputPGSEDiffusion(biExp=biExp)
        self.diffusionSummary(biExp=biExp)
        
        if self.performGradientCal == True: #save gradient calibrations to a file
            try: 
                txtgradcal='ImaxGx(A)={:.3f}, ImaxGy(A)={:.3f}, Gxcal(mT/m/A)={:.3f}, Gycal(mT/m/A)={:.3f} \n'.format(self.Axmax,self.Aymax,1000*self.GxCal, 1000*self.GyCal)
                fname=self.gradientDirectory + '\\GradientCal.txt'
                f=open(fname, 'w')
                f.write(txtgradcal)
                f.close()
                self.message('Gradient calibrations written to ' + fname +'\n')
            except:
                self.message('Gradient calibrations could not be written to ' + fname +'\n', bold=True, color='red')
    except:
      self.message('Diffusion Analysis Error')
      raise
   
  def diffusionAnalysisBiExpNWBW(self):
    '''Diffusion analysis: inputs gradient currents and 1d images to calibrate gradients, then eddy current correction data, and finaly diffusion data'''
    self.message('******Diffusion BiExp Processing******', bold=True, color='blue')
    self.ui.cbDiffusionModel.setCurrentIndex(self.ui.cbDiffusionModel.findText('Bi-exponential'))   #set diffusion model to biExponential
    self.ui.cbIntegrationRange.setCurrentIndex(self.ui.cbIntegrationRange.findText('Integrate fixed bandwidth'))    #set to integrate fixed bandwidth around major peak
    self.diffusionAnalysis(biExp=True)
#     self.setDiffusionAnalysisType()
#     try:
#         self.intbw=self.intbwWide
#         for grador in ['X','Y']:
#             if self.performGradientCal == True:
#                  self.inputGradientCurrentTraces()
#                  self.calculateImageWidths()
#                  self.inputEddyCurrentCorrection()
#             else:
#                 self.setGradOrientation(grador)
#             if self.performECC == True:
#                  self.inputEddyCurrentCorrection()
#             self.inputPGSEDiffusion(kurtosis=False, biExp=True)
#         self.diffusionSummary(biExp=True)
#         self.intbw=self.intbwWide
#         self.findIntWidth()
#         self.integrateData()
#         self.fitDiffusion(Kurtosis=False,biExp=True)
#         self.pw.exportPlot()    #puts a copy of the plot into the report
#     except:
#       self.message('Diffusion Analysis Error: BiExpNWBW')
#       raise
      
  def diffusionSpectra(self):
      '''opens PSGSEfile, baseline subtract, FFT, phasing, int widt, and plot signal vs bvalue'''
      self.message('******Diffusion Spectra Processing******', bold=True, color='blue')
      self.openDataFile(message='Open Diffusion PGSE file')
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive       
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax 
        self.message('Repeat {}, total phase adjust (deg)={:4.1f}'.format(j,self.PhaseArray[i,j]))
      self.plotData()  #plot data
      #m=4
      self.reddata=np.abs((self.data[:,:,0,0]))#define reduced data as a 2 dimensional array averaged over measurements
      #self.reddata=self.data[:,:,:,0].real.mean(axis=2)#define reduced data as a 2 dimensional array averaged over measurements
      self.reddata/=np.amax(self.reddata)     #normalize reduced data
      pc=[]     #pointcloud for T1, T2, D spectra, contains frequency, Si compnent amplitude, ADC
      for i in range (2,self.nPoints-2):
          data=(self.reddata[i-2,:]+self.reddata[i-1,:]+self.reddata[i,:]+self.reddata[i+1,:]+self.reddata[i+2,:])/5
          if np.amax(data)> self.spectraMinimum:
              params=params=DiffusionPGSE.initialize (bvalue=self.bValueArray, data=data, Kurtosis=False, biExp=False)
              pdicti=params[0] #parameter dictionary
              plist=params[1] #parameter list
              fitoutput = lmfit.minimize(DiffusionPGSE.model,pdicti,args=(self.bValueArray,data))
              pdict=fitoutput.params
              if pdict['ADC'].value<0.01 and pdict['Si'].value>0.01:
                  p=np.array([self.freq[i],pdict['Si'].value,1000*pdict['ADC'].value])
                  pc.append(p)
#               if pdict['ADC2'].value<0.01 and pdict['Si'].value>0.01:
#                   p=np.array([self.freq[i],pdict['Si2'].value,1000*pdict['ADC2'].value])
#                   pc.append(p)
      pc=np.array(pc)
      cmax=pc[:,1].sum()
      pc[:1]/=cmax  #normalize exponential amplitudes so integral = 1
      npoints=pc.shape[0]
      colorpoints=255*np.column_stack((pc[:,1],pc[:,1],pc[:,1]))#,np.full(npoints,0.5)))     #4 column color array
      color=np.column_stack((pc[:,1],pc[:,1],pc[:,1],np.full(npoints,0.5)))     #4 column color array
      #self.plotPointCloud(pc, color=color)
      self.pw=self.plotWindow(self)
#      cm = pg.colormap.get('CET-L3') # prepare a linear color map
      self.pw.win.setWindowTitle('Diffusion Spectra :' + self.fileName)
      self.pw.dplot.setLabel('bottom', "Frequency(Hz)")
      self.pw.dplot.setLabel('left', "ADC(10^-3 mm^2/s)")
      self.pw.show()
      self.pw.dplot.plot(pc[:,0],pc[:,2], pen=None, symbol='o',symbolSize=10, symbolBrush =[pg.mkBrush(v) for v in pc[:,1]])  #symbolBrush=(0,0,225)
      self.pw.dplot.plot(pc[:,0],pc[:,1], pen=None, symbol='o',symbolSize=0, symbolBrush =(50,0,225))  #symbolBrush=(0,0,225)
      self.pw.dplot.plot(self.freq,np.abs(self.data[:,0,0,0])/np.amax(self.data[:,0,0,0].real), pen=self.bpen) #plot spectra for reference
    
          
  def kurtosisAnalysis(self):
    '''Diffusion Kurtosis analysis'''
    self.message('******Diffusion Kurtosis Processing******', bold=True, color='blue')
    try:
        for g in ['X','Y']:
              self.inputGradientCurrentTraces()
              self.calculateImageWidths()
              self.inputEddyCurrentCorrection()
              self.inputPGSEDiffusion(kurtosis=True)
        self.diffusionSummary(kurtosis=True)
    except:
      self.message('Diffusion Analysis Error')
      raise   
     
  def HPAnalysis(self):
      '''opens T2CPMG file, baseline subtract, FFT, phasing, int widt, and plot signal vs echo time'''
      self.message('******Hyperpolarization Processing******', bold=True, color='blue')
      c=self.openDataFile(message='Open HP file')
      if c=='cancel':
        return
      self.subtractBaselines()
      self.FFTData()
      for j in range(self.nRepeats):    #coarse and fine phase adjustment of spectra 0
        i=0
        marg=np.argmax(np.abs(self.data[:,i,j,0]))
        dphase=np.angle(self.data[marg,i,j,0])
        self.data[:,:,j,0]= self.data[:,:,j,0]*np.exp(-1j*dphase)
        self.PhaseArray[:,j]+=dphase*180/np.pi
      self.findIntWidth() #find region around peak to use to maximize real part, will use last spectrat assuming it is fully recovered and positive 
      s='<b>Total phase adjust(deg): </b>'      
      for j in range(self.nRepeats):    #fine phase adjust   
        self.repeatIndex=j
        anglemax=self.finePhaseIJ(i,j,maxparam='nsum') #fine phase adjust, make real part most positive 
        self.data[:,:,j,0]*=np.exp(-1j*anglemax*np.pi/180)  #adjust all spectra with repect to the first spectra phase
        self.PhaseArray[:,j]+=anglemax
        s+='Repeat {}= {:4.1f}  '.format(j,self.PhaseArray[i,j]) 
      self.message(s)
      self.plotData(showCrossHairs=False)  #plot data
      self.addDataPlottoReport()
      self.integrateData()
#      self.pw.parWin.fitCPMGT2()
      self.pw.exportPlot()
                               
  def fitLorentzians(self):
      """Fits spectra with lorentzians """
      self.LorentzS0=np.zeros((self.nSpectra, self.nRepeats))
      self.Lorentztau=np.zeros((self.nSpectra, self.nRepeats))
      self.Lorentzfwhm=np.zeros((self.nSpectra, self.nRepeats))
      self.Lorentzphi=np.zeros((self.nSpectra, self.nRepeats))
      self.Lorentzf0=np.zeros((self.nSpectra, self.nRepeats))
      f=self.freq
      self.message ('<b>Lorentzian fit</b>')
      self.message('spectra#, Smax, f0(Hz), FWHM(Hz), T2*(ms), Phase(deg)')   
      for j in range (self.nRepeats): #integrate FID/spectra
        for i in range(self.nSpectra):
          datar=self.data[:,i,j,0].real
          datai=self.data[:,i,j,0].imag
          params=LorentzComplex.initialize (f=f, datar=datar, datai=datai, fitType='real')
          pdicti=params[0] #parameter dictionary
          plist=params[1] #parameter list
          fitoutput = lmfit.minimize(LorentzComplex.cLorentz,pdicti,args=(f,datar,datai))
          pdict=fitoutput.params
          s0=pdict['S0'].value
          self.LorentzS0[i,j]=s0
          f0=pdict['f0'].value
          self.Lorentzf0[i,j]=f0
          tau=pdict['tau'].value
          fwhm=1/(tau*np.pi)
          self.Lorentztau[i,j]=tau
          self.Lorentzfwhm[i,j]=fwhm
          phi=pdict['phi'].value
          self.Lorentzphi[i,j]=phi
          self.message(str(j*self.nSpectra+i) + "   {:.2e}   {:.2f}   {:.2f}   {:.2f}   {:.2f}".format(s0,f0,fwhm,1000.0*tau,phi))   
      self.message('Averages:  ' + "{:.4e}".format(np.average(self.LorentzS0)) +'  ' + "{:.2f}".format(np.average(self.Lorentzf0))  + '  ' + "{:.2f}".format(np.average(self.Lorentzfwhm)) +'  '+ "{:.2f}".format(1000.0*np.average(self.Lorentztau))+'  ' + "{:.2f}".format(np.average(self.Lorentzphi)))   
      return
         
  def plotB1z(self):
    self.B1=plotWindow(self)
    self.B1.win.setWindowTitle('B1(z) ' + str(self.repeatIndex)) 
    self.B1.show()
    self.B1.dplot.plot(self.coilz,self.b1z , pen=None, symbol='o',symbolSize=14,)
    self.B1.dplot.setLabel('bottom', "z",units='mm', **self.labelStyle) 
      
  def saveData(self, curve=1):
    self.dataFileName= QFileDialog.getSaveFileName(self,"Write data to File", "/home/file", "Ascii File freq/time real imag ... (*.dat)")
    if not self.dataFileName:
        return None
    fileName=str(self.dataFileName[0])

    f=open(fileName, 'w')
    f.write(self.fileName+'\n')
    if self.ui.cbDataType.currentText()=='FID (s)':
        xcol=self.tntfile.fid_times()
        xlabel='time(s)'
    if self.ui.cbDataType.currentText()=='Spectra (Hz)':
        xcol=self.freq
        xlabel='freq(Hz)'
    f.write(xlabel + '  real imag  real imag ...'+'\n')
    for i in range(self.nPoints):
        line="{:.4f}".format(xcol[i])
        for k in range(self.nSpectra):
            for j in range(self.nRepeats):                    
                reS="{:.4f}".format(self.data[i,k,j,0].real)
                imS="{:.4f}".format(self.data[i,k,j,0].imag)
                line+=  ', ' + reS + ', ' + imS 
        f.write(line + '\n')
    f.close()
    
  def saveMessages(self):
    self.dataFileName= QFileDialog.getSaveFileName(self,"Write message box to file", "/home/file", "Bloch Simulator Files (*.dat)")
    if not self.dataFileName:
        return None
    fileName=str(self.dataFileName)
    f=open(fileName, 'w')
    f.write(self.ui.txtMessages.toPlainText())
    f.close()

  def saveMessagestoPDF(self):
    filename = QFileDialog.getSaveFileName(self, 'Save to PDF',  "/home/file", "Bloch Simulator Files (*.pdf)")
    if filename:
        printer = QPrinter(QPrinter.HighResolution)
        printer.setPageSize(QPrinter.A4)
        printer.setColorMode(QPrinter.Color)
        printer.setOutputFormat(QPrinter.PdfFormat)
        printer.setOutputFileName(filename)
        self.ui.txtMessages.document().print_(printer) 
           
  def lackOfFitTest(self, data,modeldata, d):     
      ''' Calculates ratio of varience from model(SSLF)) to varience from noise (SSPE)
      data= 2d array of measured data with M values and d repeats, modeldata is the expected values from the model, d=# of parmeters in fit'''
      N=data.size   #total number of points= independent parameter values * repeats
      M=data.shape[0]   #independent parameter values
      Sav=np.mean(data,axis=1)[:,np.newaxis]
      SSLF=np.sum(np.square(Sav-modeldata))     #lack of fit large if the fit is bad
      SSPE=np.sum(np.square(data-Sav))          #error due to noise
 
      MSLF=SSLF/(M-d)
      MSPE=SSPE/(N-M)
      Fstat=MSLF/MSPE       #large number means model not working
      pvalue=1-fdis.cdf(Fstat,M-d,N-M)
#       print('SSLF, SSPE, N,M, dofn, dofd, Fstat,pvalue',SSLF, SSPE, N,M, M-d, N-M, Fstat, pvalue)
#       x=np.linspace(0,6,100)
#       print (fdis.cdf(x,10,10))
      return Fstat, pvalue
    
  def view2d(self):
    self.dataimage=TNMRviewer(self)
    self.dataimage.imv.setImage(np.real(self.data[:,:,:,0]))
    self.dataimage.win.show()
    self.dataimage.imv.activateWindow() 
      
  def plotBlocalPointCloud(self):
        self.updateParameters()  #update all parameters in case user modified
        self.setSpinPacketPositions() #sets spin packet position in m relative to sample tube, the sample tube position can also vary by self.sampleOffset
        bl=self.setblocal()
        c = (bl - np.min(bl))/np.ptp(bl)    #local field is encoded as 0 to 1
        color=np.column_stack((np.zeros(self.nSpinPackets),1-c,c,np.full(self.nSpinPackets,0.5)))
        self.plotPointCloud(self.SPRarray,color=color)
                
  def plotPointCloud(self, pc, color=(0,0,1,0.5)):
        self.view3Dwin = gl.GLViewWidget()
        #self.view3Dwin.opts['distance'] = 0.01
        self.view3Dwin.resize(800,800)
        self.view3Dwin.setWindowTitle('3D View ' )
        self.view3Dwin.show()
        ax = gl.GLAxisItem()
        self.view3Dwin.addItem(ax) 
        self.spPoints = gl.GLScatterPlotItem(pos=pc, color=color, size=5, pxMode=True)
        self.view3Dwin.addItem(self.spPoints)
        self.view3Dwin.show()
  
  def view3d(self):
      '''creates 3d rendering of data'''
      return
      if not hasattr(self,"view3Dwin"):   
        self.view3Dwin = gl.GLViewWidget()
        self.view3Dwin.opts['distance'] = 300
        self.view3Dwin.resize(800,800)
        self.view3Dwin.setWindowTitle('3D View ' )
      self.view3Dwin.show()
      try:
         self.view3Dwin.removeItem(self.image3DVol)
      except:
        pass
      ax = gl.GLAxisItem()
      self.view3Dwin.addItem(ax)
#       g = gl.GLGridItem()
#       g.scale(10, 10, 10)
#       self.view3Dwin.addItem(g) 
      data=np.real(self.data.astype(float))
      # /float(self.data.max())  #normalize data to 1
      d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
      d2[..., 0] = data * self.view3DColor.red()
      d2[..., 1] = data * self.view3DColor.green()
      d2[..., 2] = data * self.view3DColor.blue()
      d2[..., 3] = (data)**self.view3DTransparency * 255.   #sets transparency  
      d2[:, 0:3, 0:3] = [255,0,0,20]   #draw axes at corner of box 
      d2[0:3, :, 0:3] = [0,255,0,20]
      d2[0:3, 0:3, :] = [0,0,255,20]    
      self.image3DVol=gl.GLVolumeItem(d2)
      self.image3DVol.translate(-128,-128,-128)
      self.view3Dwin.addItem(self.image3DVol)
      #self.view3Dwin.update(self.geometry())      
      #self.view3Dwin.repaint(self.geometry())

  def openGradientCalFile(self):
    '''open gradient calibration file eg: ImaxGx(A)=7.082, ImaxGy(A)=7.232, Gxcal(mT/m/A)=48.546, Gycal(mT/m/A)=47.197 ''' 
    f = QFileDialog.getOpenFileName(self,"Open gradient calibration file", "", "Calibration Files (*.txt)")
    self.gradCalFileName=f[0]
    if  self.gradCalFileName=='':  #if cancel is pressed return
        return None
    f = open(self.gradCalFileName, 'r')
    for line in f:
        line = line.rstrip('\n; ').split('#')[0]   #remove new line and comments
        if( len( line ) >= 1):
                if line.find('ImaxGx(A)=') >= 0 :
                    self.Axmax=float(line.split('ImaxGx(A)=')[-1].split(',')[0])
                    self.message('ImaxGx(A) reset to ' + "{:.4f}".format(self.Axmax))
                if line.find('ImaxGy(A)=') >= 0 :
                    self.Aymax=float(line.split('ImaxGy(A)=')[-1].split(',')[0])
                    self.message('ImaxGy(A) reset to ' + "{:.4f}".format(self.Aymax))  
                if line.find('ImaxGz(A)=') >= 0 :
                    self.Azmax=float(line.split('ImaxGz(A)=')[-1].split(',')[0])
                    self.message('ImaxGz(A) reset to ' + "{:.4f}".format(self.Azmax))
                if line.find('Gxcal(mT/m/A)=') >= 0 :
                    self.GxCal=float(line.split('Gxcal(mT/m/A)=')[-1].split(',')[0])/1000   #/1000 because gradients are stored in T/m, gradient calibrations are in T/m/A
                    self.message('Gxcal(mT/m/A)= reset to ' + "{:.4f}".format(self.GxCal))
                if line.find('Gycal(mT/m/A)=') >= 0 :
                    self.GyCal=float(line.split('Gycal(mT/m/A)=')[-1].split(',')[0])/1000
                    self.message('Gycal(mT/m/A)= reset to ' + "{:.4f}".format(self.GyCal))  
                if line.find('Gzcal(mT/m/A)=') >= 0 :
                    self.GzCal=float(line.split('Gzcal(mT/m/A)=')[-1].split(',')[0])/1000
                    self.message('Gzcal(mT/m/A) reset to ' + "{:.4f}".format(self.GzCal))
    if self.gradOrientation=='X':   #update gradient calibrations
        self.Amax=self.Axmax
        self.GCal=self.GxCal
    if self.gradOrientation=='Y':
        self.Amax=self.Aymax
        self.GCal=self.GyCal
    if self.gradOrientation=='Z':
        self.Amax=self.Azmax
        self.GCal=self.GzCal
    self.message('<b>Set:</b> gradient direction=' + self.gradOrientation + ', Gcal(mT/m/A)=' + '{:10.3f}'.format(1000*self.GCal) +', Imax(A)=' + '{:10.3f}'.format(self.Amax))
    self.ui.leGradCal.setText('{:10.3f}'.format(1000*self.GCal))
    self.ui.leImax.setText('{:10.4f}'.format(self.Amax))
    f.close()
                      
  def message(self, m, date=False, report=True, color='black', bold=False): 
    '''prints message in gui message box'''
    '''Report=True indicates that message is to be included in output report'''
    #self.ui.txtMessages.setLineWrapMode(1)
    m= self.formatText(m, color=color, bold=bold)
    if date== True:
        self.ui.txtMessages.append(time.strftime("%c") + ': ' + m)
    else:
        self.ui.txtMessages.append(m)
    if report:    #add message to output report
      self.reportList.append(Paragraph(str(m), self.reportstyles["Normal"]))
      self.reportList.append(Spacer(1,0.1*inch)) 
    self.ui.txtMessages.verticalScrollBar().setValue(self.ui.txtMessages.verticalScrollBar().maximum())

  def formatText(self, s, color='black', bold=False):
    s='<font color=' +color + '>' + s + '</font>'
    if bold:
      s='<b>' + s + '</b>'
    return  s  #returns string with color font bold etc
    
  def addtoReport(self,m,space=0.1):
      self.reportList.append(Paragraph(str(m), self.reportstyles["Normal"]))
      self.reportList.append(Spacer(1,space*inch))
      
  def changeBackground(self):
      col = QColorDialog.getColor()
      self.dataReal.setBackground(background=col)
      self.dataImg.setBackground(background=col)

  def toggleBackground(self):
        self.plotBackgroundBlack= not self.plotBackgroundBlack
        self.setPlotBackground(black=self.plotBackgroundBlack)
        
  def setPlotBackground(self, black=True):
        if black:
          cb='k'
          cf='y'
          self.labelStyle=self.bblabelStyle
          pen=self.wpen
        else:
          cb='w'
          cf='k'
          self.labelStyle=self.wblabelStyle
          pen=self.kpen
        self.dataReal.setBackground(background=cb)
        self.dataImg.setBackground(background=cb)
#        self.dataReal.showAxis('bottom')
#        self.dataReal.showAxis('left')
        self.dataReal.setTitle( color=cf, size='18pt')
        self.dataReal.setLabel('bottom',  **self.labelStyle)
        self.dataReal.setLabel('left', **self.labelStyle)
        self.dataReal.plotItem.getAxis('left').setPen(pen)
        self.dataReal.plotItem.getAxis('bottom').setPen(pen)
        self.dataImg.setTitle( color=cf, size='18pt') 
        self.dataImg.setLabel('bottom',  **self.labelStyle)
        self.dataImg.setLabel('left',  **self.labelStyle)
        self.dataImg.plotItem.getAxis('left').setPen(pen)
        self.dataImg.plotItem.getAxis('bottom').setPen(pen)

  def calculateBaselines(self):
      '''Calculates baselines assuming end of the FID/sprectra should be 0
      caclulates from self.blRegionStart to self.blRegionStart 
      Note Tecmag arbitrarily zeros the last 1% of the data as part of their filtering'''
      bstart=int(self.blRegionStart*self.data.real.shape[0])   #calculate baselines as the average of the data beyond self.blRegion
      bstop=int(self.blRegionStop*self.data.real.shape[0])   #calculate baselines as the average of the data beyond self.blRegion
      self.dataBaseline=np.zeros((self.data.real.shape[1],self.data.real.shape[2])) + 1j*np.zeros((self.data.imag.shape[1],self.data.imag.shape[2]))    #array of complex baseline values
      #smax=np.amax(np.absolute(self.data.real))
      for i in range(self.nSpectra):
        for j in range(self.nRepeats):
          self.dataBaseline[i,j]=np.average(self.data[bstart:bstop,i,j,0].real)+1j*np.average(self.data[bstart:bstop,i,j,0].imag)     # baseline is the averagrion value of last part of the data

      
  def subtractBaselines(self, quiet=True):
      '''Calculates , subtracts baselines, then plots data'''
      self.calculateBaselines()
      for i in range(self.nSpectra):
        for j in range(self.nRepeats):
          self.data[:,i,j,0]=self.data[:,i,j,0]-self.dataBaseline[i,j]
          self.data[-self.TechMagEndZeros:,i,j,0]=0     #zero the last set of points on the  waverform since TechMag sets them to 0
      if quiet==False:
        for j in range(self.nRepeats):
          self.message('<b>Subtract baselines:</b>'+np.array2string(self.dataBaseline[:,j], precision=2, separator=',',suppress_small=True))
      self.plotData()
        
  def addPlottoReport(self,plot): 
      self.reportList.append(plot)
      self.reportList.append(Spacer(0.5,0.1*inch))
                      
  def addDataPlottoReport(self):
        '''adds both real and imaginary data plots from the main window into the output report after cahnging background to white'''
        self.setPlotBackground(black=False)
        imfile=self.reportRealImage.replace('Z',str(self.reportRealImageZ))
        exporter = pyqtgraph.exporters.ImageExporter(self.dataReal.centralWidget)
        exporter.export(imfile)  #(toBytes=True)
        im=Image(imfile, 6.5*inch, 2.0*inch)
        self.addPlottoReport(im)
        imfile=self.reportImgImage.replace('Z',str(self.reportImgImageZ))
        exporter = pyqtgraph.exporters.ImageExporter(self.dataImg.centralWidget)
        exporter.export(imfile)  #(toBytes=True)
        im=Image(imfile, 6.5*inch, 2.0*inch)
        self.addPlottoReport(im)
        self.reportRealImageZ +=1
        self.reportImgImageZ +=1
        self.setPlotBackground(black=self.plotBackgroundBlack)
        
  def printFullTNTHeader(self):
      """Expose members of the TMAG and TMG2 structures as attributes"""
      self.infowin=self.InfoWindow(self)
      self.infowin.show()
      self.infowin.info.setText(str(self.fileName))
      for name in self.tntfile.TMAG.dtype.names:
          self.infowin.info.append(name + '= '+ str(self.tntfile.TMAG[name])) 
      for name in self.tntfile.TMG2.dtype.names:
          self.infowin.info.append(name + '= '+ str(self.tntfile.TMG2[name]))
      self.infowin.info.append('TNT default tables; ' + str(self.tntfile.DELAY))
      
    #Imaging
  def displayImage(self):
      '''Makes an TNMR image window and intiates with current data'''
      self.imw=TNMRviewer(self)
      self.imw.activateWindow()
      self.imw.setWindowTitle('Image Magnitude')
      self.imw.show()
      try:
          data=np.transpose(self.data[:,:,:,:], (1, 0, 2, 3))       #Flip first and second indicies to give slice, RO, PE, Param diemntions
          self.imw.tntData=data
          self.imw.addPlotData(self.imw.tntData, imageorient=(0,1,2,3))
      except:
          raise


    
  class InfoWindow(QMainWindow):
    def __init__(self, pw, parent = None):
      '''Define info window/textbox, rv is the parent ROIView window'''    
      super(QMainWindow, self).__init__()
      self.resize(800,600)
      self.info = QTextEdit()
      self.setCentralWidget(self.info)
      self.setWindowTitle('Info')
      self.pw=pw
      
      
#   class ImageWindow(QMainWindow):
#       def __init__(self, parWin, parent = None):
#         '''Define image stack windows and menus, pw is the parent  window'''    
#         super(QMainWindow, self).__init__()
#         self.win = QMainWindow()
#         #self.win.setWindowFlags(Qt.WindowStaysOnTopHint)
#         self.win.resize(800,600)
#         self.parWin=parWin
#         self.imv = pg.ImageView( view = pg.PlotItem())
#         self.win.setCentralWidget(self.imv)
#         self.win.setWindowTitle('Image Stack')
#         point = self.rect().topRight()
#         self.win.move(point + QPoint(int(self.width()/2), 0)) 
#         self.statusBar = QStatusBar()
#         self.win.setStatusBar(self.statusBar)
#         self.menu = self.win.menuBar()
#         self.view3DColor = QColor(255, 255 ,255 , alpha=10)
#         self.view3DBackground = QColor(155, 155 ,255 , alpha=10)
#         self.view3DTransparency = 2   #set transparency scaling for 3dview, 1 = transparency set by voxel value
#         self.view3Dinvert = False    #flag to invert contrast in 3d image
#         self.InstitutionName='NIST MIG'
#         self.FoVX=160       #field of view in mm
#         self.FoVY=160
#         self.xscale=1       #=1 if image dispalys voxels, =FoV/#voxels if distance dispaly
#         self.yscale=1 
#
#         self.imageMenu = self.menu.addMenu('&Images')
#         self.actionplotMag = QAction('Plot Raw Magnitude', self)
#         self.actionplotMag.setStatusTip('Plot raw magnitude images')
#         self.imageMenu.addAction(self.actionplotMag)
#         self.actionplotMag.triggered.connect(self.plotMag)
#         self.actionplotPhase = QAction('Plot Raw Phase', self)
#         self.actionplotPhase.setStatusTip('Plot raw phase images')
#         self.imageMenu.addAction(self.actionplotPhase)
#         self.actionplotPhase.triggered.connect(self.plotPhase)
#         self.actionplotFFTMag = QAction('Plot Reconstructed Magnitude', self)
#         self.actionplotFFTMag.setStatusTip('Plot FFT magnitude images')
#         self.imageMenu.addAction(self.actionplotFFTMag)
#         self.actionplotFFTMag.triggered.connect(self.plotFFTMag)
#         self.actionplotFFTPhase = QAction('Plot Reconstructed Phase', self)
#         self.actionplotFFTPhase.setStatusTip('Plot FFT phase images')
#         self.imageMenu.addAction(self.actionplotFFTPhase)
#         self.actionplotFFTPhase.triggered.connect(self.plotFFTPhase)
#
#         self.actionWriteDICOM = QAction('Save as DICOM', self)
#         self.imageMenu.addAction(self.actionWriteDICOM)
#         self.actionWriteDICOM.triggered.connect(self.writeDicomFiles)
#
#
# #         self.actionExportImage = QAction('Export Image', self)
# #         self.imageMenu.addAction(self.actionExportImage)
# #         self.actionExportImage.triggered.connect(self.exportImages)
#
#         self.imageMenu = self.menu.addMenu('&Processing')
#         self.actionResizeImage = QAction('Resize Image', self)
#         self.imageMenu.addAction(self.actionResizeImage)
#         self.actionResizeImage.triggered.connect(self.interpolate)
#         self.actionNormalizeImage = QAction('normalize Images', self)
#         self.imageMenu.addAction(self.actionNormalizeImage)
#         self.actionNormalizeImage.triggered.connect(self.normalizeImages)
#         self.actionSetFoV = QAction('Set field of view (FoV)', self)
#         self.imageMenu.addAction(self.actionSetFoV)
#         self.actionSetFoV.triggered.connect(self.setFoV)
#         self.actionToggleScaledView = QAction('Toggle axes from Voxel index <-> distance(mm)', self)
#         self.imageMenu.addAction(self.actionToggleScaledView)
#         self.actionToggleScaledView.triggered.connect(self.toggleScaledVoxelView)
#
#         self.imageMenu = self.menu.addMenu('&3D Images')
#         self.action3DImage = QAction('Plot 3D Reconstructed Image', self)
#         self.imageMenu.addAction(self.action3DImage)
#         self.action3DImage.triggered.connect(self.view3d)
#         self.scaledData=False       #flag to indicate if scaled or raw data should be shown
#
#         self.imageMenu = self.menu.addMenu('&NIST_7T_MRI')
#         self.actionProcessScout = QAction('Process Scout', self)
#         self.imageMenu.addAction(self.actionProcessScout)
#         self.actionProcessScout.triggered.connect(self.processScout) 
#         self.actionProcessGEFlash = QAction('Process GE_Flash or SE_T1', self)
#         self.imageMenu.addAction(self.actionProcessGEFlash)
#         self.actionProcessGEFlash.triggered.connect(self.processGEFlash)               
#         self.actionExtractCh1 = QAction('Extract Ch1', self)
#         self.imageMenu.addAction(self.actionExtractCh1)
#         self.actionExtractCh1.triggered.connect(self.extractCh1)
#         self.actionExtractCh2 = QAction('Extract Ch2', self)
#         self.imageMenu.addAction(self.actionExtractCh2)
#         self.actionExtractCh2.triggered.connect(self.extractCh2)
#         self.actionExtractCh3 = QAction('Extract Ch3', self)
#         self.imageMenu.addAction(self.actionExtractCh3)
#         self.actionExtractCh3.triggered.connect(self.extractCh3)
#         self.actionExtractCh4 = QAction('Extract Ch4', self)
#         self.imageMenu.addAction(self.actionExtractCh4)
#         self.actionExtractCh4.triggered.connect(self.extractCh4)
#         self.actionSeparateImages = QAction('Separate Images', self)
#         self.imageMenu.addAction(self.actionSeparateImages)
#         self.actionSeparateImages.triggered.connect(self.separateImages)        
#
#         self.imv.getView().setLabel('bottom',"H","voxel index")
#         self.imv.getView().setLabel('left',"V","voxel index")
#         self.dataAxes={'t':0, 'x':1, 'y':2}
#
#         self.inf1 = pg.InfiniteLine(movable=True, angle=90, label='x={value:0.2f}', 
#                        labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
#         self.inf2 = pg.InfiniteLine(movable=True, angle=0, pen=(0, 0, 200),  hoverPen=(0,200,0), label='y={value:0.2f}', 
#                    labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
#         self.proxy2 = pg.SignalProxy(self.imv.view.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMove)
#
#
#       def setFoV(self):
#         '''Sets field of view and e-display ReconMag image with new FoVs'''
#         fov, OK= QInputDialog.getDouble(self,"Input field of view in X-direction", "FoVX(mm)", value=self.FoVX, decimals=1)
#         if OK:
#             self.FoVX=fov
#         fov, OK= QInputDialog.getDouble(self,"Input field of view in Y-direction", "FoVY(mm)", value=self.FoVY, decimals=1)
#         if OK:
#             self.FoVY=fov
#         self.xscale=self.FoVX/self.fftdata.shape[-2]
#         self.yscale=self.FoVY/self.fftdata.shape[-1]
#         self.imv.getView().setLabel('bottom',"X","mm")
#         self.imv.getView().setLabel('left',"Y","mm")
#         self.plotFFTMag()
#
#       def addCrossHairs(self): 
#           self.imv.addItem(self.inf1)
#           self.imv.addItem(self.inf2)
#
#       def exportImages(self):
#         self.imv.export("geeks")
#
#       def toggleScaledVoxelView(self):
#           if self.xscale==1 and self.yscale==1:
#               self.xscale=self.FoVX/self.fftdata.shape[-2]
#               self.yscale=self.FoVY/self.fftdata.shape[-1]
#               self.imv.getView().setLabel('bottom',"X","mm")
#               self.imv.getView().setLabel('left',"Y","mm")
#           else:
#               self.xscale=1
#               self.yscale=1
#               self.imv.getView().setLabel('bottom',"H","voxel index")
#               self.imv.getView().setLabel('left',"V","voxel index")
#           self.plotFFTMag()
#
#       def plotMag(self):
#         self.win.setWindowTitle('Raw Image Magnitude, shape={}, dataMax={:.2f}, fftMax={:.2f}'.format(self.data.shape,np.max(self.data), np.max(self.fftdata)))
#         self.imv.getView().setLabel('bottom',"H","voxel index")
#         self.imv.getView().setLabel('left',"V","voxel index")
#         self.imv.setImage(np.absolute(self.data),axes=self.dataAxes)
#         self.rawDataType='RawMag'
#       def plotPhase(self):
#         self.win.setWindowTitle('Raw Image Phase, shape={}, dataMax={:.2f}, fftMax={:.2f}'.format(self.data.shape,np.max(self.data), np.max(self.fftdata)))
#         self.imv.getView().setLabel('bottom',"H","voxel index")
#         self.imv.getView().setLabel('left',"V","voxel index")
#         self.imv.setImage(np.angle(self.data),axes=self.dataAxes)
#         self.rawDataType='RawPhase'
#       def plotFFTMag(self):
#         self.win.setWindowTitle('Reconstructed Image Magnitude, shape={}, dataMax={:.2f}, fftMax={:.2f}'.format(self.data.shape,np.max(self.data), np.max(self.fftdata)))
#         self.imv.setImage(np.absolute(self.fftdata),axes=self.dataAxes,scale = (self.xscale,self.yscale))
#         self.rawDataType='ReconMag'
#       def plotFFTPhase(self):
#         self.win.setWindowTitle('Reconstructed Image Phase, shape={}, dataMax={:.2f}, fftMax={:.2f}'.format(self.data.shape,np.max(self.data), np.max(self.fftdata)))
#         self.imv.setImage(np.angle(self.fftdata),axes=self.dataAxes,scale = (self.xscale,self.yscale))
#         self.rawDataType='ReconPhase'
#
#       def addPlotData(self, data):
#             '''expect data to be 2d image or a stack of 2d images with first dimenson being the index'''
#             self.data=data
#             dt=np.fft.fftshift(data,axes=(-2,-1))
#             self.fftdata=np.fft.fftshift(np.fft.fft2(dt,axes=(-2, -1)),axes=(-2,-1))
#             dshape=self.data.shape
#             dmax=np.max(self.data)
#             fftmax=np.max(self.fftdata)
#             self.win.setWindowTitle('Reconstructed Image Magnitude, shape={}, dataMax={:.2f}, fftMax={:.2f}'.format(dshape,dmax, fftmax))
#             self.imv.setImage(np.absolute(self.fftdata),axes=self.dataAxes)
#             self.rawDataType='ReconMag'
#
#       def gaussianFilter(self, data):
#             self.data=data
#             self.fftdata=np.fft.fftshift(np.fft.fft2(data),axes=(1,2))
#             self.win.setWindowTitle('Reconstructed Image Magnitude')
#             self.imv.setImage(np.absolute(self.fftdata),axes=self.dataAxes)
#
#       def view3d(self):
#         '''creates 3d rendering of data'''  
#         self.view3Dwin = gl.GLViewWidget()
#         self.view3Dwin.opts['distance'] = 300
#         self.view3Dwin.resize(800,800)
#         self.view3Dwin.setWindowTitle('3D View ' )
#         self.view3Dwin.show()
#         try:
#           self.view3Dwin.removeItem(self.image3DVol)
#         except:
#           pass
#         ax = gl.GLAxisItem()
#         self.view3Dwin.addItem(ax)
#         if self.scaledData:
#             data=np.absolute(self.scaledImage.astype(float))
#         else:
#             data=np.absolute(self.fftdata.astype(float))
#         data=10*data/np.amax(data)
#         d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
#         d2[..., 0] = data * self.view3DColor.red()
#         d2[..., 1] = data * self.view3DColor.green()
#         d2[..., 2] = data * self.view3DColor.blue()
#         d2[..., 3] = (data)**self.view3DTransparency * 255.   #sets transparency  
#         d2[:, 0:3, 0:3] = [255,0,0,20]   #draw axes at corner of box 
#         d2[0:3, :, 0:3] = [0,255,0,20]
#         d2[0:3, 0:3, :] = [0,0,255,20]    
#         self.image3DVol=gl.GLVolumeItem(d2)
#         #self.image3DVol.translate(-128,-128,-128)
#         self.view3Dwin.addItem(self.image3DVol)
#
#       def interpolate(self):
#         '''scales and interpolates 2d images using PIL.image.resize, changes image array size
#         3d image are scaled using scipy.ndimage.zoom'''
#         interp= ['bicubic', 'nearest', 'bilinear', 'lanczos']
#         scalex,ok = QInputDialog.getDouble(self, "Scale X",'the number of rows will be increased (or decreased) by', value=2.0, min=0.1, max=10.0,decimals=2)
#         if not ok:
#           return
#         scaley,ok = QInputDialog.getDouble(self, "Scale Y",'the number of columns will be increased (or decreased) by', 2.0, 0.1, 10.0,1)
#         if not ok:
#           return
#         scalez,ok = QInputDialog.getDouble(self, "Scale Z",'the number of images in stack will be increased (or decreased) by', 2.0, 0.1, 10.0,1)
#         if not ok:
#             return
#         self.scaledImage=zoom(np.absolute(self.fftdata), (scalez,scalex,scaley),order=3)
#         self.scaledPhase=zoom(np.angle(self.fftdata), (scalez,scalex,scaley),order=3)
# #        self.pw.message('Scaled image stack using scipy.ndimage.zoom, scalex={:.2f}, scaley={:.2f}, scalez={:.2f} ,interpolation=cubic spline image shape= {},{},{} \n'.format(scalex,scaley,scalez,s[0],s[1],s[2]))
#         self.win.setWindowTitle('Reconstructed Image Magnitude, scaled by sx={}, sy={}, sz={}'.format(scalex, scaley,scalez))
#         self.imv.setImage(np.absolute(self.scaledImage),axes=self.dataAxes)
#         self.scaledData=True
#
#       def normalizeImages(self):
#           max=np.amax(self.data, axis=(1,2))
#
#       def extractCh1(self):
#           self.addPlotData(self.data[:,0::4,:])
#       def extractCh2(self):
#           self.addPlotData(self.data[:,1::4,:])
#       def extractCh3(self):
#           self.addPlotData(self.data[:,2::4,:])
#       def extractCh4(self):
#           self.addPlotData(self.data[:,3::4,:])
#
#       def separateImages(self,n=3):
#           '''Separates 3 consecutive images in echo, used for tri-scouts'''
#           if n==0:
#               n=3
#           dshape=self.data.shape
#           columns=int(dshape[1]/n)
#           rows=dshape[2]
#           self.data.shape=(3, columns, rows)
#           self.addPlotData(self.data)
#
#       def processScout(self): 
#           self.data=self.data[:,0::4,:]     #extract channel 1 from  CH1,2,3,4 data 
#           dshape=self.data.shape
#           columns=int(dshape[1]/3)  #Saggital, coronal, axial need to be separated 
#           rows=dshape[2]
#           self.data.shape=(3, columns, rows)
#           self.parWin.message('<b>msMRI Scout processing:</b> Extracted CH1, separated Coronal,Saggital,Axial slices')
#           if rows< columns:
#               npad=int((columns-rows)/2)
#               self.data=np.pad(self.data, ((0,0), (0,0), (npad, npad)), 'constant', constant_values=((0, 0),(0, 0),(0, 0)))
#               self.parWin.message('Zero padded phase encode by {}'.format(2*npad))
#           self.addPlotData(self.data)
#
#       def processGEFlash(self): 
#           self.data=self.data[:,0::4,:]     #extract channel 1 from  CH1,2,3,4 data 
#           dshape=self.data.shape 
#           columns=dshape[1]
#           rows=dshape[2]
#           self.parWin.message('<b>msMRI GE processing:</b> Extracted CH1')
#           if rows< columns:
#               npad=int((columns-rows)/2)
#               self.data=np.pad(self.data, ((0,0), (0,0), (npad, npad)), 'constant', constant_values=((0, 0),(0, 0),(0, 0)))
#               self.parWin.message('Zero padded phase encode by {}'.format(2*npad))
#           self.addPlotData(self.data)
#
#       def writeDicomFiles(self, filename):
#             nimages=self.data.shape[0]
#             ncolumns=self.data.shape[1]
#             nrows=self.data.shape[2]
#             f = QFileDialog.getSaveFileName(self,'Enter DICOM filename', '', "DICOM Files (*.dcm)")
#             if self.rawDataType=='ReconMag':
#                     dData=np.absolute(self.fftdata)
#             if self.rawDataType=='ReconPhase':
#                     dData=np.angle(self.fftdata)   
#             if self.rawDataType=='RawMag':
#                     dData=np.absolute(self.data)
#             if self.rawDataType=='RawPhase':
#                     dData=np.angle(self.data)   
#             if f[0]=='':
#                 return 'cancel'
#             self.DICOMfileName=str(f[0])      #make sure fileName is a string, not a Qstring        
#             for i in range(nimages):     #write out images in separate dicom files
#                 pixel_array=dData[i,:,:]
#                 fileName = self.DICOMfileName.replace('.dcm', str(i) + ".dcm")
#                 fileName = fileName.replace('/','\\')        
#                 # Populate required values for file meta information
#                 file_meta = Dataset()
#                 file_meta.MediaStorageSOPClassUID = b'Secondary Capture Image Storage'
#                 file_meta.MediaStorageSOPInstanceUID = b'1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
#                 file_meta.ImplementationClassUID = b'1.3.6.1.4.1.9590.100.1.0.100.4.0'
#                 # Create the FileDataset instance (initially no data elements, but file_meta supplied)
#                 ds = FileDataset(fileName, {}, file_meta=file_meta, preamble=b"\0"*128)
#                 ds.Modality = b'MR'
#                 ds.ContentDate = str(datetime.date.today()).replace('-','')
# #                 ds.ContentTime = str(time.time()) #milliseconds since the epoch
#                 ds.StudyInstanceUID =    b'1.3.6.1.4.1.9590.100.1.1.124313977412360175234271287472804872093'
#                 ds.SeriesInstanceUID = b'1.3.6.1.4.1.9590.100.1.1.369231118011061003403421859172643143649'
#                 ds.SOPInstanceUID =        b'1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
#                 ds.SOPClassUID = b'MR Image Storage'
#                 ds.SecondaryCaptureDeviceManufctur = b'Python 2.7.3'
#             ## These are the necessary imaging components of the FileDataset object.
#                 ds.SamplesPerPixel = 1
#                 ds.PhotometricInterpretation = b"MONOCHROME2"
#                 ds.PixelRepresentation = 0
#                 ds.HighBit = 15
#                 ds.BitsStored = 16
#                 ds.BitsAllocated = 16
#                 ds.SmallestImagePixelValue = b'\\x00\\x00'
#                 ds.LargestImagePixelValue = b'\\xff\\xff'
#                 #Add MRI data elements 
# #                ds.bValue= self.bValue[i]
#                 ds.Columns=ncolumns
#                 ds.Rows= nrows                
#                 ds.PixelSpacing=[self.FoVX/ncolumns,self.FoVY/nrows]
# #                ds.FoVX=self.FoVX
# #                ds.FoVY=self.FoVY
# #                 ds.ImageOrientationPatient=self.ImagePosition[i].tolist()
#                 ds.InstitutionName=self.InstitutionName
# #                 ds.PixelBandwidth= self.PixelBandwidth[i] 
# #                 ds.PixelSpacing =[self.PixelSpacingY[i],self.PixelSpacingX[i]]         
#
# #                 ds.PatientName = self.PatientName[i]
# #                 ds.PatientID = "123456"
# #                 ds.ProtocolName=self.ProtocolName[i]
# #                 ds.RepetitionTime = self.TR[i]
# #                 ds.SeriesDescription=self.SeriesDescription[i]
# #                 ds.EchoTime = self.TE[i]
# #                 ds.FlipAngle= self.FA[i]
# #                 ds.InversionTime=self.TI[i] 
# #                 ds.ImageOrientationPatient[3:6]=self.ColumnDirection[i].tolist() 
# #                 ds.ImageOrientationPatient[:3]=self.RowDirection[i].tolist()
# #                 ds.MagneticFieldStrength = self.MagneticFieldStrength[i]
# #                 ds.Manufacturer = self.Manufacturer[i]            
# #                 pixel_array=np.transpose(self.PA[i])        #numpy    pixel array
#                 if self.rawDataType== "RawPhase" or self.rawDataType== "ReconPhase":
#                     pixel_array=(pixel_array+np.pi) *10000    # phase data -pi to pi is converted to (0 to 2pi)*1000
#                     pixel_array = pixel_array.astype(np.uint16)
#                     #print "Adjusting phase to 16 bit integer"
#                 if self.rawDataType== "RawMag" or self.rawDataType== "ReconMag":
#                     scale=65536/np.amax(pixel_array)
#                     pixel_array*=scale       #rescale so that max value is 16 bits                             
#                     pixel_array = pixel_array.astype(np.uint16)
#                 ds.PixelData = pixel_array.tobytes()    # image byte data
#                 # Set the transfer syntax
#                 ds.is_little_endian = True
#                 ds.is_implicit_VR = True
#                 pydicom.filewriter.write_file(fileName, ds, write_like_original=False)
#
#       def mouseMove(self,evt): 
#             '''mouse move event to display location and values'''
#             pos = evt[0]  ## using signal proxy turns original arguments into a tuple
#             try:
#                 if self.imv.view.sceneBoundingRect().contains(pos):
#                     mousePoint = self.imv.view.vb.mapSceneToView(pos)
#                     self.Xindex=int(mousePoint.x()) 
#                     self.Yindex=int(mousePoint.y() )
#                     im=self.imv.getImageItem().image
#                     if self.Xindex>=0 and self.Xindex<im.shape[0] and self.Yindex>=0 and self.Yindex<im.shape[1]:
#                         self.statusBar.showMessage('X={}, Y={}, Value={:.2f}'.format(self.Xindex, self.Yindex, im[self.Xindex,self.Yindex]),5000)#
#             except:
#                 pass   

  class generalPlot(QMainWindow):
      def __init__(self,parWin, image=None, parent=None):
            '''Defines window for plotting noise data'''    
            super(QMainWindow, self).__init__()
            self.win=self   #QMainWindow()
            self.dplot = pg.PlotWidget()
            self.dplot.enableAutoRange(x=True,y=True)
            self.plotName=''
            self.win.setCentralWidget(self.dplot)
            self.win.resize(1200,750)
            self.win.setWindowTitle('Acquired Data')
            self.parWin=parWin  #parent window or calling form
            
            self.menu = self.win.menuBar() 
            
  class plotNoise(QMainWindow):
      def __init__(self,parWin, image=None, parent=None):
            '''Defines window for plotting noise data'''    
            super(QMainWindow, self).__init__()
            self.win=self   #QMainWindow()
            self.dplot = pg.PlotWidget()
            self.dplot.enableAutoRange(x=True,y=True)
            self.plotName=''
            self.win.setCentralWidget(self.dplot)
            self.win.resize(1200,750)
            self.win.setWindowTitle('Acquired Data')
            self.parWin=parWin  #parent window or calling form
            self.ypen=pg.mkPen('y', width=2)
            
            self.menu = self.win.menuBar() 
            
            self.imageMenu = self.menu.addMenu('&Analysis')
            self.actionTimeMax = QAction('Set maximum time', self.win)
            self.imageMenu.addAction(self.actionTimeMax)
            self.actionTimeMax.triggered.connect(self.setMaxTime)
            self.actionAveragingParameters = QAction('Set Averaging parameters', self.win)
            self.imageMenu.addAction(self.actionAveragingParameters)
            self.actionAveragingParameters.triggered.connect(self.changeSavGolParameters)

            self.imageMenu = self.menu.addMenu('&Plot')
            
            self.actionPlotPhase = QAction('Plot Phase', self.win)
            self.imageMenu.addAction(self.actionPlotPhase)
            self.actionPlotPhase.triggered.connect(self.plotPhase)
            
            self.actionPlotField = QAction('Plot Field Noise', self.win)
            self.imageMenu.addAction(self.actionPlotField)
            self.actionPlotField.triggered.connect(self.plotField)
                       
            self.actionPlotFieldSpectra = QAction('Plot Field Spectra Data', self.win)
            self.imageMenu.addAction(self.actionPlotFieldSpectra)
            self.actionPlotFieldSpectra.triggered.connect(self.plotFieldSpectra)
                
 
            self.tMax=1
            self.freq=0
            self.aveWindow=10
            self.savgolPolynomialOrder=3
            self.fmin=-400      #plot freeqeuncy min in Hz
            self.fmax=400
            
      def plotFieldSpectra(self, noiseav=8):
          #self.ftNoiseData=fft.fftshift(fft.fft(np.unwrap(self.phaseData[:self.maxIndex,0,0])))
          self.ftNoiseData=fft.fftshift(fft.fft(self.field))/self.field.shape[0]
          self.dplot.clear()
          self.dplot.setLabel('bottom',"Frequency(Hz)")
          self.dplot.setLabel('left',"Field Spectra (nT)")
          self.dplot.setTitle("Field Spectra",size='14pt', color=(100,200,200))
          fmax=0.5/self.dwellTime
          npoints=len(self.field)
          df=fmax/npoints
          freqs=2*np.arange(npoints)*df-fmax
          #fieldspectra=gaussian_filter1d(np.absolute(self.ftNoiseData)*freqs, 8)
          self.dplot.plot(freqs,np.absolute(self.ftNoiseData)*1E9, name='Noise Spectra', pen=self.ypen) 
          self.dplot.setXRange(self.fmin, self.fmax, padding=0)
          
      def plotPhase(self):
          self.dplot.clear()
          self.dplot.setLabel('bottom',"Time(s)")
          self.dplot.setLabel('left',"Phase (rad)")
          self.dplot.setTitle("Phase",size='14pt', color=(100,200,200))
          self.phase=np.unwrap(self.phaseData[:self.maxIndex])
          self.field=(self.maxIndex/self.tMax)*savgol_filter(self.phase, self.aveWindow,self.savgolPolynomialOrder, deriv=1)/self.parWin.Gamma
          self.dplot.plot(self.tfid[:self.maxIndex],self.phase, name='Phase', pen=self.ypen) 
          
      def plotField(self):
          '''Plots magnetic field derived from taking the smoothed derivative of the FID phase and dividing by gyromagnetic ratio'''
          self.dplot.clear()
          self.dplot.setLabel('bottom',"Time(s)")
          self.dplot.setLabel('left',"Field (nT)")
          self.dplot.setTitle("Field",size='14pt', color=(100,200,200))
          
          self.dplot.plot(self.tfid[:self.maxIndex],self.field*1E9, name='Field Noise', pen=self.ypen) 
                    
      def setMaxTime(self, tMax=0):
          '''Gets upper index corresponding to tMax of desired phase plot, data above this index is ignored'''
          if tMax==0:
              tMax, OK= QInputDialog.getDouble(self,"Max time", "Input maximum time (s)", 1, 0, 100)
              if not OK:
                  return
          self.tMax=tMax
          difference_array = np.absolute(self.tfid-tMax)
          self.maxIndex = difference_array.argmin()
          self.plotPhase()
              
      def changeSavGolParameters(self):
          avwindow, OK= QInputDialog.getInt(self,"Averaging window size", "Input window size", self.aveWindow, 2, 100)
          if OK:
              self.aveWindow=avwindow
          polyorder, OK= QInputDialog.getInt(self,"Averaging polynomioal ordere", "Input polynomial order", self.savgolPolynomialOrder, 1, 6)
          if OK:
              self.savgolPolynomialOrder=polyorder
              self.plotPhase()
                                      
  class plotWindow(QMainWindow):
      def __init__(self,parWin, image=None, parent=None):
            '''Defines window for plotting data, fit to data, residuals. Stores lists of data, fits, and residuls'''    
            super(QMainWindow, self).__init__()
            self.win=self   #QMainWindow()
            self.dplot = pg.PlotWidget()
            self.plotName=''
            self.win.setCentralWidget(self.dplot)
            self.win.resize(1200,750)
            self.win.setWindowTitle('Acquired Data')
            self.parWin=parWin  #parent window or calling form
            self.menu = self.win.menuBar()
            self.penr = QPen(Qt.red, 0.1)
            self.peny = QPen(Qt.yellow, 0.1)
            self.penb = QPen(Qt.blue, 0.1)
            self.kpen=pg.mkPen('k', width=3)
            self.bblabelStyle = {'color':'w', 'font-size': '18px'}
            self.bbtitleStyle = {'color':'w', 'font-size': '18px'}
            self.wblabelStyle = {'color':'k', 'font-size': '18px'}
            self.wbtitleStyle = {'color':'k', 'font-size': '18px'}
            self.labelStyle=self.bblabelStyle
            self.titleStyle=self.bbtitleStyle
            self.symbolSize=12
            self.logX=False
            self.logY=False
            self.plotBackgroundBlack=True
            self.setPlotBackground(self.plotBackgroundBlack)
            self.bPlotAll=True       #flag to determine whether to plot all curves or just a selected one.
            self.bClearData=True     #flag to determine whether to clear old data before plotting new data
            self.plotTitle=''
            self.dplot.showAxis('right')
            self.dplot.showAxis('top')
#            self.dplot.setLabel('bottom',"Time(s),**self.labelStyle)
#            self.dplot.setLabel('left',"Signal","AU",**self.labelStyle)
            self.symb=['o', 's', 'd', 't', 't1', 't2','t3', 'p','+', 'h', 'star','x']
            self.xData=[]
            self.yData=[]       #primary y data
            self.xResData=[]       #list of residuals
            self.yResData=[]
            self.xFitData=[]
            self.yFitData=[]
            self.header=[]
            self.Resheader=[]
            self.Fitheader=[]
            self.dataDescr=[]
            self.fitStats=[]
            self.xSNR=[]
            self.ySNR=[]
            self.dataType=''    #type of NMR data FID, diffusion, nutation, T1IR
            self.rawDataType='' #type of raw data real imadinairy or magnitude pahse for plots, NMR data lways stored as complex values
            self.dataStart=0        # start fits at integer dataStart
            self.dataStop=-1         # stop fits at dataStop 
            
            self.imageMenu = self.menu.addMenu('&File')    
            self.actionSaveData = QAction('Save data', self.win)
            self.imageMenu.addAction(self.actionSaveData)
            self.actionSaveData.triggered.connect(self.saveData)
            self.actionReadData = QAction('Read data', self.win)
            self.imageMenu.addAction(self.actionReadData)
            self.actionReadData.triggered.connect(self.readData)
            
            self.imageMenu = self.menu.addMenu('&Analysis')    
            self.actionT1IR = QAction('T1 IR fit', self.win)
            self.imageMenu.addAction(self.actionT1IR)
            self.actionT1IR.triggered.connect(self.T1Fit)
            self.actionCPMGT2 = QAction('CPMG/SE T2 fit', self.win)
            self.imageMenu.addAction(self.actionCPMGT2)
            self.actionCPMGT2.triggered.connect(parWin.fitCPMGT2)
            self.actionCPMGT2biEx = QAction('CPMG/SE bi-exp T2 fit', self.win)
            self.imageMenu.addAction(self.actionCPMGT2biEx)
            self.actionCPMGT2biEx.triggered.connect(parWin.fitCPMGT2biEx)
            self.actionT1rho = QAction('T1 rho fit', self.win)
            self.imageMenu.addAction(self.actionT1rho)
            self.actionT1rho.triggered.connect(parWin.fitT1rho)
            
            self.actionDiffusion = QAction('Diffusion fit', self.win)
            self.imageMenu.addAction(self.actionDiffusion)
            self.actionDiffusion.triggered.connect(self.fitDiffusion)
            self.actionDiffusionTruncated = QAction('Truncated Diffusion fit', self.win)
            self.imageMenu.addAction(self.actionDiffusionTruncated)
            self.actionDiffusionTruncated.triggered.connect(self.fitDiffusionTruncated)
            self.actionDiffusionBiExp = QAction('Diffusion BiExp fit', self.win)
            self.imageMenu.addAction(self.actionDiffusionBiExp)
            self.actionDiffusionBiExp.triggered.connect(self.fitDiffusionBiExp)
            self.actionDiffusionKurtosis = QAction('Diffusion/ Kurtosis fit', self.win)
            self.imageMenu.addAction(self.actionDiffusionKurtosis)
            self.actionDiffusionKurtosis.triggered.connect(self.fitDiffusionKurtosis)
            self.actionNutation = QAction('Nutation fit', self.win)
            self.imageMenu.addAction(self.actionNutation)
            self.actionNutation.triggered.connect(parWin.fitNutation)
            self.actionFitAll = QAction('Fit separately', self.win)
            self.imageMenu.addAction(self.actionFitAll)
            self.actionFitAll.triggered.connect(self.fitAll)
            self.actionRefData = QAction('Save as refeence data', self.win)
            self.imageMenu.addAction(self.actionRefData)
            self.actionDiffusion.triggered.connect(self.refData)
            self.actionCalData = QAction('Divide by reference data', self.win)
            self.imageMenu.addAction(self.actionCalData)
            self.actionDiffusion.triggered.connect(self.calData)
            
            self.imageMenu = self.menu.addMenu('&Plot')
            self.actionPlotData = QAction('Plot Data', self.win)
            self.imageMenu.addAction(self.actionPlotData)
            self.actionPlotData.triggered.connect(self.plotData)    
            self.actionPlotResiduals = QAction('Plot Residuals', self.win)
            self.imageMenu.addAction(self.actionPlotResiduals)
            self.actionPlotResiduals.triggered.connect(self.plotResiduals)
            self.actionPlotSNR = QAction('Plot SNR', self.win)
            self.imageMenu.addAction(self.actionPlotSNR)
            self.actionPlotSNR.triggered.connect(self.plotSNR)

            self.imageMenu = self.menu.addMenu('&Data')
            self.actionScale0to1 = QAction('Scale data 0 to 1', self.win)
            self.imageMenu.addAction(self.actionScale0to1)
            self.actionScale0to1.triggered.connect(self.scaleData0to1)    

                        
            self.imageMenu = self.menu.addMenu('&Plot Options')
            self.actionPlotAll = QAction('Plot Single Curve', self.win)
            self.imageMenu.addAction(self.actionPlotAll)
            self.actionPlotAll.triggered.connect(self.togglePlotAll)
            self.actionClearPlot = QAction('Clear Plot', self.win)
            self.imageMenu.addAction(self.actionClearPlot)
            self.actionClearPlot.triggered.connect(self.cleardata)       
            self.actionClearData = QAction('Clear Data Off', self.win)
            self.imageMenu.addAction(self.actionClearData)
            self.actionClearData.triggered.connect(self.toggleClearData)    
            self.actionToggleBackground = QAction('Toggle Background', self.win)
            self.imageMenu.addAction(self.actionToggleBackground)
            self.actionToggleBackground.triggered.connect(self.toggleBackground)
            self.actionChangeBackground = QAction('Change Background', self.win)
            self.imageMenu.addAction(self.actionChangeBackground)
            self.actionChangeBackground.triggered.connect(self.changeBackground)
            self.actionAddTitle = QAction('Add to Title', self.win)
            self.imageMenu.addAction(self.actionAddTitle)
            self.actionAddTitle.triggered.connect(self.addToTitle)
            self.actionAddCursor = QAction('Add Cursor', self.win)
            self.imageMenu.addAction(self.actionAddCursor)
            self.actionAddCursor.triggered.connect(self.addCursor)
            self.actionSetXlog = QAction('Toggle X-Axis Log', self.win)
            self.imageMenu.addAction(self.actionSetXlog)
            self.actionSetXlog.triggered.connect(self.setXLogMode)
            self.actionSetYlog = QAction('Toggle Y-Axis Log', self.win)
            self.imageMenu.addAction(self.actionSetYlog)
            self.actionSetYlog.triggered.connect(self.setYLogMode)
            self.actionSetSymbolSize = QAction('Set Symbol Size', self.win)
            self.imageMenu.addAction(self.actionSetSymbolSize)
            self.actionSetSymbolSize.triggered.connect(self.setSymbolSize)            
                        
            self.imageMenu = self.menu.addMenu('&Report')
            self.actionAddPlottoReport = QAction('Add Plot to Report', self.win)
            self.imageMenu.addAction(self.actionAddPlottoReport)
            self.actionAddPlottoReport.triggered.connect(self.exportPlot)
            #try: 
            self.inf1 = pg.InfiniteLine(movable=True, angle=90, label='x={value:0.4e}', 
                       labelOpts={'position':0.1, 'color': (200,200,100), 'fill': (200,200,200,50), 'movable': True})
            self.inf2 = pg.InfiniteLine(movable=True, angle=0, pen=(0, 0, 200),  hoverPen=(0,200,0), label='y={value:0.4e}', 
                       labelOpts={'color': (200,0,0), 'movable': True, 'fill': (0, 0, 200, 100)})
            #except:
            #    pass

      def setPlotTitle(self):
          self.dplot.setTitle(self.plotTitle,size='14pt', color=(100,200,200))
          
      def setSymbolSize(self):
            s, OK= QInputDialog.getInt(self,"Input symbol size", "", min=0, max=40)
            if OK:
                self.symbolSize= s
      
      def setLogMode(self, logx=False, logy=False):
          self.logX=logx
          self.logY=logy
          self.dplot.setLogMode(self.logX, self.logY)
          
      def setXLogMode(self):
          self.logX= not self.logX
          self.setLogMode(logx=self.logX,logy=self.logY)

      def setYLogMode(self):
          self.logY= not self.logY
          self.setLogMode(logx=self.logX,logy=self.logY)
                     
      def togglePlotAll(self):
            '''Toggles between plotting all curves and just selected curve'''
            self.bPlotAll=not self.bPlotAll
            if self.bPlotAll:
              self.actionPlotAll.setText('Plot Single Curve')
            else:
              self.actionPlotAll.setText('Plot All Curves')
              
      def toggleClearData(self):
          self.bClearData=not self.bClearData
          if self.bClearData:
              self.actionClearData.setText('Clear Data Off')
          else:
              self.actionClearData.setText('Clear Data On')     
      
      def dataType(self,datatype):
          self.dataType=datatype
                     
      def refData(self):
          '''save a data a a reference'''
          pass 
      
      def calData(self):
          '''save a data a a reference'''
          pass 
            
      def fitDiffusion(self):
          self.parWin.fitDiffusion(dataType=self.dataType)
          
      def fitDiffusionTruncated(self):
          '''Fits truncated data and plots data and fit'''
          self.dataStart+=1
          self.cleardata()
          self.plotData()
          self.parWin.message('Truncated diffusion data start= {}, stop={}'.format(self.dataStart, self.dataStop), color='orange')
          self.parWin.fitDiffusion(dataType=self.dataType, dataStart=self.dataStart)
                    
      def fitDiffusionBiExp(self):
          self.parWin.fitDiffusion(dataType=self.dataType, biExp=True) 
                   
      def fitDiffusionKurtosis(self):          
          self.parWin.fitDiffusion(dataType=self.dataType, Kurtosis=True)
                    
      def T1Fit(self):
            self.parWin.dataInt=np.zeros((len(self.xData[0]),len(self.yData)))            
            self.parWin.TI=self.xData[0]
            for j in range (self.parWin.nRepeats): #
                self.parWin.dataInt[:,j]=self.yData[j]
                #self.parWin.dataInt[:,j]=self.yData[j]
                #self.parWin.dataInt[:,j]=self.yData[j]
            self.parWin.fitT1IR()
               
      def readData(self):
        self.dataFileName= QFileDialog.getOpenFileName(self,"Read data file", "/home/file", "Bloch Simulator Files (*.dat *.csv)")
        if not self.dataFileName:
            return None
        data=np.transpose(np.genfromtxt(str(self.dataFileName), delimiter=',', comments = '#'))
        self.xData.append(data[0][np.isfinite(data[0])])
        self.yData.append(data[1][np.isfinite(data[1])])
        self.xData.append(self.xData[0])
        self.yData.append(data[2][np.isfinite(data[2])])
        self.xData.append(self.xData[0])
        self.yData.append(data[3][np.isfinite(data[3])])
        self.plotData()
      
      def saveData(self):
        self.dataFileName= QFileDialog.getSaveFileName(self,"Write data to File", "/home/file", "Bloch Simulator Files (*.dat)")
        if not self.dataFileName:
            return None
        header=str('pyNMR data, window title= ' + self.win.windowTitle())
        nel=max(len(self.xData[0]),len(self.xFitData[0]))
        s='time(s) \n' 
        for i  in range(nel):
            if i < len(self.xData[0]):
                s += str(self.xData[0][i]) + ','
            else:
                s += 'time, '       
            for d in self.yData: 
                if i < len(d):
                    s +=   str(d[i]) +', '
                else:
                    s += 'signal, '
            for j,d in enumerate(self.xResData):
                if j==0: #only print first x column, all should be the same 
                    if i < len(d):
                        s +=  "{:.5f}".format(d[i])+ ', ' 
                    else:
                        s += 'tres, '
            for d in self.yResData: 
                if i < len(d):
                    s +=  "{:.5f}".format(d[i]) + ', ' 
                else:
                    s += 'res, '
            for j,d in enumerate(self.xFitData): 
                if j==0: #only print first x column, all should be the same
                    if i < len(d):
                        s +=  str(d[i]) +', ' 
                    else:
                        s += 'tFit, ' 
            for d in self.yFitData: 
                if i < len(d):
                    s +=  str(d[i])+ ', ' 
                else:
                    s += 'yFit, ' 
            s += '\n'                                  
        f= open(self.dataFileName, 'w')
        f.write(header)
        f.write(s)
        f.close()

    
      def addData(self,x,y, description='', header='', title=''):
            self.xData.append(x)
            self.yData.append(y)
            self.dataDescr.append(description)
            self.header.append(header)
            
      def addResData(self,x,y, description='', header='', title=''):
            '''add residuals to the data set'''
            self.xResData.append(x)
            self.yResData.append(y)
            #self.dataDescr.append(description)
            self.Resheader.append(header) 
            
      def addFitData(self,x,y, description='', header='', title=''):
            self.xFitData.append(x)
            self.yFitData.append(y)
            self.dataDescr.append(description)
            self.Fitheader.append(header)
            
      def addSNRData(self,x,y, description='', header='', title=''):
            self.xSNR.append(x)
            self.ySNR.append(y)
            self.dataDescr.append(description)

                                               
      def addFitStats(self,f):
          self.fitStats.append(f)
          
      def cleardata(self):
          self.dplot.clear()
          
      def fitAll(self):
          self.parWin.fitAll=not self.parWin.fitAll
          if self.parWin.fitAll:
              QAction.setText('"Fit separately')
          else:
              QAction.setText('"Fit all curves')
              
      def changeBackground(self):
            col = QColorDialog.getColor()
            self.dplot.setBackground(background=col) 
            
      def plotSNR(self):
          self.dplot.clear()
          self.parWin.dplot.setLogMode(False, False)
          for i in range(len(self.xSNR)):
              self.dplot.plot(self.xSNR[i], self.ySNR[i], pen=None, symbol=self.symb[i],symbolSize=self.symbolSize)
        
      def plotResiduals(self):
          if self.bClearData:
              self.dplot.clear()
          self.dplot.setLogMode(False, False)
          if self.bPlotAll:
              for i in range(len(self.xResData)):
                p=pg.mkPen(pg.intColor(i), width=2)
                self.dplot.plot(self.xResData[i], self.yResData[i], pen=p, symbol=self.symb[i%12],symbolSize=self.symbolSize, symbolPen=p)
          if not self.bPlotAll:
                pn, OK= QInputDialog.getInt(self,"Input curve number", "", min=0, max=len(self.xResData)-1)
                if OK:
                    p=pg.mkPen(pg.intColor(pn), width=2)
                    self.dplot.plot(self.xResData[pn], self.yResData[pn], pen=p, symbol=self.symb[pn%12],symbolSize=self.symbolSize,symbolPen=p)
                       
      def plotData(self):
          if self.bClearData:
              self.dplot.clear()
          if self.bPlotAll:
              for i in range(len(self.xData)):
                  self.dplot.plot(self.xData[i][self.dataStart:self.dataStop], self.yData[i][self.dataStart:self.dataStop], pen=None, symbol=self.symb[i%12],symbolSize=self.symbolSize)
          if not self.bPlotAll and len(self.xData)>0:
                    pn, OK= QInputDialog.getInt(self,"Input curve number", "", min=0, max=len(self.xData)-1)
                    if OK:
                        p=pg.mkPen(pg.intColor(pn), width=2)
                        self.dplot.plot(self.xData[pn][self.dataStart:self.dataStop], self.yData[pn][self.dataStart:self.dataStop], pen=p, symbol=self.symb[pn%12],symbolPen=p, symbolSize=self.symbolSize) 
#           for i in range(len(self.xFitData)):
#               self.dplot.plot(self.xFitData[i], self.yFitData[i], pen=self.penb, symbol=None,symbolSize=24,) 

      def resetDataRange(self):
          self.dataStart=0
          self.dataStop=-1
          
      def scaleData0to1(self):
          '''linearly scale y data so first point is 1 and last point is zero'''
          self.dplot.clear()
          self.dplot.setLogMode(False, True)
          if self.bPlotAll:
              for i in range(len(self.xData)):
                  self.dplot.plot(self.xData[i], (self.yData[i]-self.yData[i][-1])/(self.yData[i][0]-self.yData[i][-1]), pen=None, symbol=self.symb[i%12],symbolSize=self.symbolSize)
              for i in range(len(self.xFitData)):
                  self.dplot.plot(self.xFitData[i], (self.yFitData[i]-self.yFitData[i][-1])/(self.yFitData[i][0]-self.yFitData[i][-1]), pen=self.peny)

                                  
      def addToTitle(self): 
          s, OK= QInputDialog.getText(self,"Append to Title", "")
          if OK:
              #l=pg.TextItem(s,color=(255,0,0),anchor=(0.5,0.5))
              #l.setPos(10,10)
              #self.dplot.addItem(l)
              self.plotTitle += '\n'  + s
              self.setPlotTitle()
      
      def addCursor(self): 
          self.dplot.addItem(self.inf1)
          self.dplot.addItem(self.inf2)
              
      def exportPlot(self):
          '''Sends plot to pdf output report'''
          self.setPlotBackground(black=False)
          app.processEvents()   #seems to be needed otherwise axes screwed up
          imfile=self.parWin.reportFitImage.replace('Z',str(self.parWin.reportFitImageZ))   #name image file with integer indicate plot number in the report
          exporter = pyqtgraph.exporters.ImageExporter(self.dplot.centralWidget)  #self.dplot.centralWidget
          exporter.parameters()
          exporter.export(imfile)  #(toBytes=True) save image to jpg file
          im=Image(imfile, 5.5*inch, 3.5*inch)  #open jpg file
          self.parWin.addPlottoReport(im)
          self.parWin.reportFitImageZ+=1 
          self.setPlotBackground(black=self.plotBackgroundBlack)
          
      def toggleBackground(self):
        self.plotBackgroundBlack= not self.plotBackgroundBlack
        self.setPlotBackground(black=self.plotBackgroundBlack)
        
      def setPlotBackground(self,black=True):
        '''changes the plot background/title and label colors from dark/light to light/dark'''
        if black:
          cb='k'
          cf='w'
          self.labelStyle=self.bblabelStyle
          pen=self.parWin.wpen
        else:
          cb='w'
          cf='k'
          self.labelStyle=self.wblabelStyle
          pen=self.kpen
        xlabel=self.dplot.getAxis('bottom').labelText
        ylabel=self.dplot.getAxis('left').labelText
        self.dplot.setBackground(background=cb)
        self.dplot.setLabel('bottom',xlabel,**self.labelStyle)
        self.dplot.setLabel('left',ylabel, **self.labelStyle)
        self.dplot.plotItem.getAxis('left').setPen(pen)
        self.dplot.plotItem.getAxis('bottom').setPen(pen)
 
class textFormat:
   '''unicode text formatting: textFormat.Purple' will turn text purple'''
   Purple = '\033[95m'
   Cyan = '\033[96m'
   DarkCyan = '\033[36m'
   BlueE = '\033[94m'
   Green = '\033[92m'
   Yellow = '\033[93m'
   Red = '\033[91m'
   Bold = '\033[1m'
   Underline = '\033[4m'
   End = '\033[0m'                          

#useful code templates
#p, ok =  QInputDialog.getDouble(self, 'Phase adjust', 'Enter phase(deg)')
#Useful for debugging Qt applications where the app closes without giving error message
sys._excepthook = sys.excepthook 
def exception_hook(exctype, value, traceback):
    print("Missed Exception:", exctype, value, traceback)
    sys._excepthook(exctype, value, traceback) 
    sys.exit(1) 
#*******
sys.excepthook = exception_hook 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    clipboard=app.clipboard()
    app.setStyleSheet("QWidget{font-size: 8pt;}") 
    main = pyNMR()
    main.show()
    sys.exit(app.exec_())
