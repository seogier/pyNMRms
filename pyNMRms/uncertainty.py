# -*- coding: utf-8 -*-
"""
Created on Jan 27, 2018
NMR Uncertainty calculation
Uses uncertaintyGui created from uncertaintyGui.ui by QT4
  execute from system shell to regenerate GUIs
    designer\\pyuic4 designer\\uncertaintyGUI.ui -o pyNMRms\\uncertaintyGUI4.py
    designer\\pyuic5 designer\\uncertaintyGUI.ui -o pyNMRms\\uncertaintyGUI5.py 

author: serus
"""
import sys
import numpy as np
from pyqt import *         #imports required PyQt modules, tries PyQT4 then PyQt5
if pyqtVersion==4:
    from uncertaintyGUI4 import Ui_Uncertainty
if pyqtVersion==5:
    from uncertaintyGUI5 import Ui_Uncertainty
from scipy import constants


class UncertaintyWindow(QMainWindow):
  def __init__(self ,pw,parent = None):
        super(UncertaintyWindow, self).__init__()
        self.ui = Ui_Uncertainty()
        self.ui.setupUi(self)
        self.pw=pw    #sets access to parent window name space
        self.setWindowTitle('Uncertainty')
        self.ui.pbCalculateUncertainty.clicked.connect(self.calculateUncertainty)
  
  def calculateUncertainty(self):
        self.pw.T1Uncertainty=0.015*self.pw.T1reported
        self.pw.T2Uncertainty=0.007*self.pw.T2reported
        self.ui.leT1Uncertainty.setText(str(self.pw.T1Uncertainty))
        self.ui.leT1Uncertainty.setText(str(self.pw.T1Uncertainty))
        self.pw.message(self.pw.formatText('*'*100,color='blue'))
        t1='T1(ms)=' + "{:.3f}".format(self.pw.T1reported*1000) + ', T1 uncertainty(ms)=' + "{:.3f}".format(self.pw.T1Uncertainty*1000)
        t1+=', @ B(T)=' + "{:.3f}".format(self.pw.magneticField) + ', Temperature(C)=' + "{:.3f}".format(self.pw.temperature)
        self.pw.message(self.pw.formatText(t1, color='blue', bold=True))
        t2='T2(ms)=' + "{:.3f}".format(self.pw.T2reported*1000) + ', T2 uncertainty(ms)=' + "{:.3f}".format(self.pw.T2Uncertainty*1000)
        t2+=', @ B(T)=' + "{:.3f}".format(self.pw.magneticField) + ', Temperature(C)=' + "{:.3f}".format(self.pw.temperature) + ', TauCPMG(ms)=' + "{:.3f}".format(self.pw.tauCPMG * 1000)
        self.pw.message(self.pw.formatText(t2, color='blue', bold=True))
        qc=self.pw.formatText('All QC tests=', color='blue', bold=True)
        if self.pw.allQCtests=='Pass':
          c='green'
        else:
          c='red'
        qc+=self.pw.formatText(self.pw.allQCtests, color=c, bold=True)
        self.pw.message(qc)  
        self.pw.message(self.pw.formatText('*'*100,color='blue'))       
        
        
if __name__ == '__main__':
    pass