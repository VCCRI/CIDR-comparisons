## Make the Simulation folder the working directory.
from datetime import datetime
from ZIFA import ZIFA
from ZIFA import block_ZIFA
import numpy
Y=numpy.loadtxt("results/forZifa.csv",delimiter=",",skiprows=1)
startTime = datetime.now()
Z2, model_params = block_ZIFA.fitModel(Y, 2)
print datetime.now() - startTime
numpy.savetxt("results/Z2.csv",Z2,delimiter=",")
X=Y

