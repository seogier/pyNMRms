'''
Each model is referred to using a modelname and must contain must contain three methods
  intialize
  model
  fit
DiffusionPGSE : Simple pulsed gradient spin echo exponential decay model to extract diffusion coeficient 
bvalue in s/mm^2, typically 0 to 10000 s/mm^2
ADC in mm^2/s
last modification: 9-10-17
'''

import lmfit
import numpy as np


def initialize (bvalue=None, data=None, B=0,Kurtosis=False, biExp=False):    
    """initialize parameters for Diffusion model
    Can fit simple exponential, biexponential or kurtosis models"""
    nparams =3      #max number of parameters, some may be fixed

    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    guess=1/np.mean(bvalue)
    params.add('ADC', value= guess, min=0, vary = True)
    paramlist.append('ADC')
    params.add('K', value= 0, min=-10, vary = Kurtosis)
    paramlist.append('K')
    params.add('Si', value= np.amax(data), vary = True)
    paramlist.append('Si')
    params.add('ADC2', value= guess/5, min=0, vary = biExp)
    paramlist.append('ADC2')
    if biExp:
        params.add('Si2', value= np.amax(data)/10, vary = True)
    else:
        params.add('Si2', value= 0, vary = False)
    paramlist.append('Si2')
    params.add('B',  value= 0,  min=0., vary = False)
    paramlist.append('B')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def model(params, bvalue, data):
    """ Diffusion PGSE model; """
    B = params['B'].value
    Si = params['Si'].value
    ADC = params['ADC'].value
    Si2 = params['Si2'].value
    ADC2 = params['ADC2'].value
    K = params['K'].value
    x=bvalue*ADC
    x2=bvalue*ADC2
    model = Si*np.exp(-x+K*x**2/6)+ Si2*np.exp(-x2)+B
    return (model - data)

def fit(params, bvalue, data):
    """fits signal vs bvalue data to diffusion model to extract apparent diffusion coefficient ADC and Kurtosis K"""
    result = lmfit.minimize(model, params, args=(bvalue, data))
    final = data + result.residual
    return final

