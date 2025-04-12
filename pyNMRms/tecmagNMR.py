'''
Created on Sep 4, 2017

@author: serus
'''

class tecmagNMR(object):
    '''
    Parameters describing NIST tecmagNMR
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        self.birdCageWidth=14.0       #birdcageRFcoil WIDTH
        self.B1amplitude=4.0E-6            #typical RF amplitude in T