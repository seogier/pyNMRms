# pyNMRms
 
## Requirements
Anaconda3       #this includes python,  scipy, numpy, and lots of other useful packages;  download and run the anaconda installer for the system that you are using
pyqtgraph            this includes all of the plotting packages used, use the latest version, current version is 
lmfit                this includes the nonlinear least squares fitting
PyOpenGL             used for 3d plotting, not critical
nmrglue              opens Varian and Bruker data
reportlab            generates pdf reports
pydicom==1.4.2       for importing and exporitng dicom files, not critical, only package at present that needs old version
 
 Can run the following lines in anaconda3\scripts where pip is located   
    pip install pyqtgraph           
    pip install lmfit               
    pip install PyOpenGL            
    pip install nmrglue             
    pip install reportlab            
    pip install pydicom==1.4.2        
     
Current package versions
Python=3.12.4 | packaged by Anaconda, Inc. | (main, Jun 18 2024, 15:03:56) [MSC v.1929 64 bit (AMD64)];
PyQt=5.15.10;
pyqtgraph=0.13.7
## Execution

From the `pyNMRms\` directory, run `pyNMR.py`.