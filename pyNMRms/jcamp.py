'''
Created on Apr 22, 2022

@author: russek
'''

class MyClass(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
     
     def convertNMR_RawPicoSpin80_to_ndarray(txt):
#where txt is the raw text of the full .jdx file.
tfile = tempfile.mkstemp(prefix="nmr_", suffix=".jdx")[1]
ffile = open(tfile, 'w')
ffile.write(txt)
ffile.close()
dct, raw = ng.jcampdx.read(tfile)
data = numpy.empty(len(raw[0]), dtype='complex128')
data.real = raw[0]; data.imag = raw[1]
#process fid (time-domain) to spectrum (frequency-domain)
data = ng.proc_base.zf_size(data, int(dct['$PSZEROFILLING'][0])) #8192
data = ng.proc_base.fft(data)
data = ng.proc_base.ps(data, p0=float(dct['$PSPHASECORRECTION'][0].replace(' degrees', '')), p1=0.0)
#data = ng.proc_base.di(data)
data = numpy.absolute(data)
#determine the ppm scale
udct = ng.jcampdx.guess_udic(dct, data)
udct[0]['car'] = 0.
udct[0]['sw'] = float(dct['$PSMAXFREQ.TOPLOT'][0].replace(' Hz', '')) - float(dct['$PSMINFREQ.TOPLOT'][0].replace(' Hz', '')) #4000. #3985.
uc = ng.fileiobase.uc_from_udic(udct)
return uc.ppm_scale(), data   