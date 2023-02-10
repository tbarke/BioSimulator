import numpy as np
import cell
import enviornment
import utils
import simulation
import configuration
import MICalc
import numpy
import random
import randomClass
import matplotlib.pyplot as plt
import math
import pickle
#import numpy as np
import gc
import run
import time
import runSuite
from scipy import special
from scipy.stats import binom
from matplotlib.colors import LogNorm

c = configuration.configuration()
c.runOutputFlags.save = True
c.simParams.simLength = 10
c.runOutputFlags.compressSave = True
run.testRun(c)
exit(-1)