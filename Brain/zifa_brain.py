## Make the Brain folder the working directory.
from datetime import datetime
from ZIFA import ZIFA
from ZIFA import block_ZIFA
import numpy
Y=numpy.loadtxt("Results/forZifa.csv",delimiter=",",skiprows=1)
startTime = datetime.now()
Z4, model_params = block_ZIFA.fitModel(Y, 4)
print datetime.now() - startTime
numpy.savetxt("Results/Z4.csv",Z4,delimiter=",")
X=Y

