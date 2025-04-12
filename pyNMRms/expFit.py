'''
Each model is referred to using a modelname and must contain must contain three methods
  intialize
  model
  fit
Simple exponential fit
last modification: 9-10-17
'''

import lmfit
import numpy as np


def initialize (t=None, data=None, B=0):    
    """initialize parameters for Diffusion model"""
    nparams =3      #max number of parameters, some may be fixed

    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    guess=1/np.mean(t)
    params.add('tau', value= guess, min=1E-5, vary = True)
    paramlist.append('tau')
    params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    params.add('B',  value= 0,  min=0., vary = False)
    paramlist.append('B')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def ObjFnc(params, t, data):
    """ exp model; """
    B = params['B'].value
    Si = params['Si'].value
    tau = params['tau'].value

    model = Si*np.exp(-t/tau)+ B
    return (model - data)

def fit(params, t, data):
    """fits signal to exponential with time constant tau"""
    result = lmfit.minimize(ObjFnc, params, args=(t, data))
    final = data + result.residual
    return final

