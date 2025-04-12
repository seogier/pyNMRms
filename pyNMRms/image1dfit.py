'''
Each model is referred to using a modelname and must contain must contain three methods
  intialize
  ObjFnc
  fit
image1d : fits a 1 d image of a cylinder of diameter dc and height hc, ax specifies the orientation of the 1d slice, 0=perp to cyl axis, 1=parallel to axis
last modification: 9-10-17
'''

import lmfit
import numpy as np


def initialize (f=None, data=None, ax=1):    
    """initialize parameters"""
    nparams =5      #max number of parameters, some may be fixed

    params = lmfit.Parameters()   #define parameter dictionary
    paramlist = []    # list of parameters used for this model
    Aguess=np.max(data)     #initial guess is average of input time values
    params.add('A', value= Aguess, min=0, vary = True)     #Signal maximum
    paramlist.append('A')
    params.add('f0', value= 0, vary = True) #center frequency, hopefully around 0Hz
    paramlist.append('f0')
    params.add('gf', value= abs(np.max(f)/5), min= 0.0,vary = True)    #gradient in Hz gf=G*gamma*dc
    paramlist.append('gf')
    params.add('bl', value= 0, vary = True)     #baseline
    paramlist.append('bl')
    params.add('ax', value= ax, vary = False)     #cylinder axis
    paramlist.append('ax')
    return [params,paramlist]

# define objective function: returns the array to be minimized
def ObjFnc(params, f, data):
    """ model is signal as a function of frequency"""
    A = params['A'].value
    f0 = params['f0'].value
    gf = params['gf'].value
    bl = params['bl'].value
    ax = params['ax'].value

    if ax==0:       # If image axis is perpendicular to cylinder axis, signal is constant inside +-gf
        #s=np.piecewise(f, [f < f0-gf,np.absolute(f) <= f0+gf, f > f0+gf], [0,1,0])
        s=1-(f-f0)**8/gf**8
        s=np.clip(s,0,None)     #set negative values to 0
        #s=np.sign(s)
        model = A*s+bl
    else:   #image plane is along cylinder axis 
        s=1-(f-f0)**2/gf**2
        s=np.clip(s,0,None)     #set negative values to 0
        model = A*s**0.5+bl
    return (model - data)

def fit(params, f, data):
    """fits signal vs frequency"""
    result = lmfit.minimize(ObjFnc, params, args=(f, data))
    final = data + result.residual
    return final

