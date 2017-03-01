from datetime import datetime
from ZIFA import ZIFA
from ZIFA import block_ZIFA
import numpy
Y=numpy.loadtxt("/Users/paulyLin/Dropbox/CIDR_paper/manuscript/Codes/Pancreatic/Results/forZifa.csv",delimiter=",",skiprows=1)
startTime = datetime.now()
Z4, model_params = block_ZIFA.fitModel(Y, 4)
print datetime.now() - startTime
numpy.savetxt("/Users/paulyLin/Dropbox/CIDR_paper/manuscript/Codes/Pancreatic/Results/Z4.csv",Z4,delimiter=",")
X=Y

