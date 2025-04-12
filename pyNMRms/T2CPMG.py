'''
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
T2CPMG : Simple T2 CPMG spin echo exponential decay model 
last modification: 9-10-17
'''

import lmfit
import numpy as np


def initializeT2CPMG (TE=None, data=None, B=0):    
    """initialize parameters for T2CPMG model"""
    nT2CPMGparams =3      #max number of parameters, some may be fixed

    T2params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    T2guess=np.mean(TE)
    T2params.add('T2', value= T2guess, min=0., vary = True)
    paramlist.append('T2')
    T2params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    T2params.add('B',  value= 0,  min=0., vary = False)
    paramlist.append('B')
    return [T2params,paramlist]

# define objective function: returns the array to be minimized
def T2CPMG(params, TE, data):
    """ T2-CPMG model; TE echo time array"""
    B = params['B'].value
    Si = params['Si'].value
    T2 = params['T2'].value

    model = Si*np.exp(-TE/T2)+ B
    return (model - data)

def fitT2CPMG(params, TE, data):
    """fits signal vs TE data to T2SE model"""
    result = lmfit.minimize(T2CPMG, params, args=(TE, data))
    final = data + result.residual
    return final

