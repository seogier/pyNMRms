"""
Created on Fri Oct 11 16:30:54 2013
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T1IR : T1 inversion recovery standard model  Signal = A*(1-B * np.exp(-TI/T1))
last modification: 6-3-14
"""

import lmfit
import numpy as np


def initializeT1IR (nParam=None,TI=None, data=None):    
    """initialize parameters for T1IR absolute value model"""
    nT1IRparams =3      #max number of parameters, some may be fixed
    if nParam!=None:
        return nT1IRparams
    T1params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    T1guess=TI[np.argmin(np.absolute(data))]/np.log(2) #minimum signal should occur at ln(2)T1
    T1params.add('T1', value= T1guess, min=0, vary = True)
    paramlist.append('T1')
    T1params.add('A', value= np.amax(data), vary = True)
    paramlist.append('A')
    T1params.add('delta',  value= 1,  min=0, max=1.5, vary = True)
    paramlist.append('delta')
    return [T1params,paramlist]

# define objective function: returns the array to be minimized
def T1IR(params, TI, data):
    """ T1-IR model abs(exponential); TI inversion time array, T1 recovery time"""
    delta = params['delta'].value
    A = params['A'].value
    T1 = params['T1'].value

    model = A*(1-(1+delta) * np.exp(-TI/T1))
    return (model - data)

def fitT1IR(params, TI, data):
    """fits signal vs TI data to T1IRabs model"""
    result = lmfit.minimize(T1IR, params, args=(TI, data))
    final = data + result.residual
    return final
