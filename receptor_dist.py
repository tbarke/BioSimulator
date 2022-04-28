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
# -*- coding: utf-8 -*-

SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA (negative if random), locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25, 25]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change)]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]

CellMetaStats = [2.0, 0.0, 5, 1.0, "non", "non"]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []
CellLocationStats = [Celllocations]

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die]
EnviornmentParams = [100, False, False]

CellMetaStats[5] = "non"
simLength = EnviornmentParams[0]
simRepeat = 2
fianlTimeDivx = range(simLength)
finalTimeDiv = [0]*simLength
survivaly = []

ABreceptors, concs = run.moveTrackRunReceptor(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True, True)



#sorted_zip = sorted(zipped)
Areceptors = [0.0] * 100
Breceptors = [0.0] * 100
index = [0] * 100
for i in range(len(ABreceptors[2])):
    for j in range(100):
        if ABreceptors[2][i] == j:
            index[j] += 1
            Areceptors[j] = Areceptors[j] + (ABreceptors[0][i]-Areceptors[j])/index[j]
            Breceptors[j] = Breceptors[j] + (ABreceptors[1][i]-Breceptors[j])/index[j]

x = range(100)
print(Areceptors)

for i in range(len(Breceptors)):
    #bx[i] = bx[i] * -1
    Breceptors[i] = Breceptors[i] * -1



#fig = plt.figure()
#gs = fig.add_gridspec(2, hspace=0)
#axs = gs.subplots(sharex=True)

fig, axs = plt.subplots()
ax1 = axs.bar(range(100), Areceptors)
ax2 = axs.bar(range(100), Breceptors, color = 'm')
#axs[0].set_yticks([0.0,0.5,1.0])
plt.xticks(rotation=90)
plt.tight_layout()
# Use absolute value for y-ticks
ticks = axs.get_yticks()
axs.set_yticklabels([int(abs(tick)) for tick in ticks])
axs.set_xlabel("location")
axs.set_ylabel("receptor count")
axs.legend((ax1, ax2), ('\'A\' receptors', '\'B\' receptors'), loc='lower right')
axs.set_title("Average Receptor Count per Location (Maximum Sensitivity Cell Type)")
axs.set_xlim(0, 100)
#axs[1].set_yticks([0,150,300])
axs.grid()
plt.show()