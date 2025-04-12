'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T2CPMGbiEx : Simple T2 CPMG spin echo biexponential decay model 
last modification: 9-10-17
'''

import lmfit
import numpy as np


def initialize (TE=None, data=None, B=0):    
    """initialize parameters for T2CPMG model"""
    nparams =5      #max number of parameters, some may be fixed

    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    T2guess=np.mean(TE)     #initial guess is average of input time values
    minex=TE[0]       #set minimum time constant to 1/3 of initial time value to prevent spuri0us fast decaying components
    params.add('T2a', value= T2guess/4, min=minex, vary = True)     #try to have T2a be the short T2
    paramlist.append('T2a')
    params.add('T2b', value= T2guess, min=minex, vary = True)
    paramlist.append('T2b')
    params.add('Si', value= np.amax(data), min= 0.0,vary = True)
    paramlist.append('Si')
    params.add('Sib', value= np.amax(data),min=0.0, vary = True)
    paramlist.append('Sib')
    params.add('B',  value= 0,  min=0., vary = False)
    paramlist.append('B')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def ObjFnc(params, TE, data):
    """ T2-CPMG bi-exponential model; TE echo time array"""
    B = params['B'].value
    Si = params['Si'].value
    Sib = params['Sib'].value
    T2a = params['T2a'].value
    T2b = params['T2b'].value

    model = Si*np.exp(-TE/T2a) + Sib*np.exp(-TE/T2b)+ B
    return (model - data)

def fit(params, TE, data):
    """fits signal vs TE data to T2SE model"""
    result = lmfit.minimize(ObjFnc, params, args=(TE, data))
    final = data + result.residual
    return final

