"""
Created on 5-14-2018
Each model is referred to using a modelname and must contain must contain three methods
  intializemodelname
  modelname
  fitmodelname
CLorentz  Signal = A*(1-B * np.exp(-TI/T1))
last modification: 6-3-14
"""

import lmfit
import numpy as np

fitType='real'

def initialize (nParam=None,f=None, datar=None, datai=None, fitType='real'):    
    """initialize parameters for complex Lorentzian model
    if the fitType=real then the modelr-datar is minimized,
    if the fitType=imag then the modeli-datai is minimized"""
    nparams =4      #max number of parameters, some may be fixed
    if nParam!=None:
        return nparams
    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    S0guess=datar[np.argmax(np.absolute(datar))] #S0 is the peak height of the real part, assume phi~0
    f0guess=f[np.argmax(np.absolute(datar))]      #guess the center of lorentzian is at real data maximum
    params.add('S0', value= S0guess, vary = True)
    paramlist.append('S0')
    params.add('f0', value= f0guess, vary = True)   #f0 is the position of the lorentzian
    paramlist.append('f0')
    params.add('tau',  value= 0.1,  min=0.0, max=10.0, vary = True) #tau is T2* or FWHM=2/tau
    paramlist.append('tau')
    params.add('phi',  value= 0.0,  min=-90.0, max=90.0, vary = True) #phi is the phase in degrees
    paramlist.append('phi')
    fitType=fitType
    return [params,paramlist]

# define objective function: returns the array to be minimized
def cLorentz(params, f, datar, datai):
    """ complex Lorentzian"""
    S0 = params['S0'].value
    f0 = params['f0'].value
    tau = params['tau'].value
    tauinv=1/tau
#    fwhm=1/(tau*np.pi)
    phi=params['phi'].value*np.pi/180.0
    df=f-f0
    den=1+(2*np.pi*df*tau)**2
    real=S0*(np.cos(phi)+2*np.pi*df*tau*np.sin(phi))/den
    imag=S0*(np.sin(phi)-2*np.pi*df*tau*np.cos(phi))/den
    if fitType=='real':
      model = real
      data=datar
    if fitType=='imag':
      model = imag
      data=datai
    return (model - data)

def fit(params, f, data):
    '''fits Lorentzian'''
    result = lmfit.minimize(cLorentz, params, args=(f, data))
    final = data + result.residual
    return final
