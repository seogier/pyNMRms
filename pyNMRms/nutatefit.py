"""
Created on Fri Oct 11 16:30:54 2013
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
Nutate: nutation simple model using a trapezoidal B1 distribution
last modification: 9-4-17
"""

import lmfit
import numpy as np
from scipy import constants
from scipy import integrate
from tecmagNMR import tecmagNMR as tNMR

gamma=constants.physical_constants["proton gyromag. ratio"][0]*(1-constants.physical_constants["proton mag. shielding correction"][0])                             
birdCageWidth=14.0

def initializeNutate (nParam=None,tau=None, data=None, B1max=100.0E-6):    
    """initialize parameters for Nutate absolute value model"""
    nNutateparams =3      #max number of parameters, some may be fixed
    if nParam!=None:
        return nNutateparams
    Nutateparams = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
                                        
    S0guess=np.maximum(data) #Set initial scale at maximum of input data
    Nutateparams.add('S0', value= S0guess, min=0, vary = True)
    paramlist.append('S0')
    Nutateparams.add('B1max', value= tNMR.B1amplitude, vary = True)
    paramlist.append('B1max')
    Nutateparams.add('dbc',  value= 5,  min=0, max=tNMR.birdCageWidth/2, vary = True)
    paramlist.append('B')
    return [Nutateparams,paramlist]

# define objective function: returns the array to be minimized
def Nutate(params, tau, data):
    """ nutate model ; tau RF pulse duration, T1 recovery time"""
    S0 = params['S0'].value
    B1max = params['B1max'].value
    dbc = params['dbc'].value

    model = signal(b1,B1max,tau)
    return (model - data)

def b1(z,b1max, d):
    w=birdCageWidth/2
    if np.absolute(z) > w +d/2:
        return 0.0
    if np.absolute(z) < w-d/2:
        return b1max
    if np.absolute(z-w) < d/2:
        return b1max*((w-z)/d+0.5)
    if np.absolute(z+w) < d/2:
        return b1max*((w+z)/d+0.5)

        
def signal(b1,b1max,tau):
    s=np.sin(gamma*b1*tau*b1max)*b1
    signal=np.trapz(s)
    return signal

def fitNutate(params, TI, data):
    """fits signal vs TI data to Nutateabs model"""
    result = lmfit.minimize(Nutate, params, args=(TI, data))
    final = data + result.residual
    return final
