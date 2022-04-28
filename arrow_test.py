from scipy import special
import cell
import enviornment
import utils
import simulation
import MICalc
import numpy
import random
import randomClass
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import math
import numpy as np
import gc
import run


def VonMises(var, maximum, offset, length):
    length = length
    arr = []
    kappa = 1 / var
    mu = 2 * math.pi / length * offset
    bessel = (special.iv(0, kappa))
    sum = 0.0
    increment = 2 * math.pi / length
    for i in range(length):
        index = i - length / 2
        x = index * increment
        numerator = math.exp(kappa * math.cos(x - mu))
        arr.append(numerator / (2 * math.pi * bessel))
        sum += numerator / (2 * math.pi * bessel)
    arr = (arr / sum) * maximum
    return arr

ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/100), 500, -25.5, 25.5, 0, 50]

print(utils.calculateBoundKinetic(20, 0.05, 2, 1))
exit(-1)

ConcentrationA = VonMises(ConcParams[7], ConcParams[8], ConcParams[9], 100)
ConcentrationB = VonMises(ConcParams[7], ConcParams[8], ConcParams[10], 100)
ConcentrationC = VonMises(ConcParams[7], ConcParams[8], ConcParams[11], 100)
x = range(100)
plt.bar(x,ConcentrationA, alpha = 0.5, label = 'A', color = 'red', edgecolor = 'red')
plt.bar(x, ConcentrationB, alpha = 0.5, label = 'B', color = 'blue', edgecolor = 'blue')
plt.bar(x,ConcentrationC, alpha = 0.5, label = 'C', color = 'green', edgecolor = 'green')
plt.legend()
plt.xlabel("Location")
plt.ylabel("Concentration")
plt.grid()
plt.show()