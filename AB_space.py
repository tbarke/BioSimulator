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
import math
import numpy as np
import gc
import run


SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA (negative if random), locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25, 25]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change)]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]

CellMetaStats = [2.0, 0.0, 9, 1.0, "non", "non"]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []
CellLocationStats = [Celllocations]

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die]
EnviornmentParams = [1000, False, False]

CellMetaStats[5] = "measured"
simLength = EnviornmentParams[0]
simRepeat = 1
fianlTimeDivx = range(simLength)
finalTimeDiv = [0]*simLength
survivaly = []


arr = run.moveTrackRunReceptorIntAB(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True, False)
Aint = arr[1][2]
Bint = arr[1][3]
divloc = arr[1][4]

print(Aint)
print(Bint)
print(divloc)

m = max(Aint)
n = max(Bint)
m = max([m,n])

xx = []
yx = []
for i in range(len(divloc)):
    if divloc[i] == 1:
        xx.append(Aint[i])
        yx.append(Bint[i])

point1 = [0, m]
point2 = [0, m]

title = "Measured"
#title = "Maximum Sensitivity"

#plt.plot(point1, point2, "black", alpha=0.6)
plt.plot(Aint, Bint, 'b' )
plt.scatter(xx, yx)
plt.grid()
plt.xlim(-2,300)
plt.ylim(-2,300)
plt.plot([0, 0], [m, m])
plt.hlines(5*5, 5*5, m, "red", alpha=0.6)
plt.vlines(5*5, 5*5, m, "red", alpha=0.6)
#plt.title("Cell Trajectory in A and B Internal Molecule\n for " + title + " Cell Type")
plt.xlabel("'A' molecules")
plt.ylabel("'B' molecules")

plt.show()
