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

scale = 1
SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA (negative if random), locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25, 25, 0, 50]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change)]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]
cell_stress_array = [1/4,1/3,1/2, 2/3, 3/4,1/2]
CellMetaStats = [1/cell_stress_array[5], 0.0, 5, 1.0, "non", "non", 0, "simulate", scale]
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

divide = []
locs  = []
divide2 = []
locs2  = []
splits= []
splits2 = []
for i in range(simRepeat):
    print(str(i) + ": Repeat Number")
    timeRun = run.moveTrackRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True, False)
    if i == 0:
        locs.append(timeRun[0])
        divide.append(timeRun[1])
        splits.append(timeRun[2])
    if i == 1:
        locs2.append(timeRun[0])
        divide2.append(timeRun[1])
        splits2.append(timeRun[2])

print(locs)
print(divide)
print(locs2)
print(splits2)
print(divide2)

x = range(100)

#print(locs)
#print(divide)

#locs = [[51, 50, 33, 24, 13, 5, 88, 76, 74, 81, 76, 93, 95, 94, 98, 9, 10, 14, 23, 27, 18, 12, 93, 84, 78, 70, 58, 57, 38, 24, 27, 26, 31, 19, 20, 9, 94, 99, 74, 64, 84, 78, 73, 71, 58, 46, 48, 40, 40, 58, 55, 51, 46, 63, 66, 60, 50, 39, 37, 41, 48, 66, 72, 63, 52, 53, 36, 22, 18, 19, 7, 82, 83, 71, 73, 78, 74, 80, 84, 79, 72, 73, 73, 80, 84, 73, 60, 65, 55, 49, 48, 18, 23, 14, 22, 30, 35, 27, 31, 37]]
#divide = [[0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]]
#locs2 = [[51, 50, 33, 24, 13, 5, 88, 76, 74, 81, 76, 93, 95, 94, 98, 9, 10, 14, 23, 27, 18, 12, 93, 84, 78, 70, 58, 57, 38, 24, 27, 26, 31, 19, 20, 9, 94, 99, 74, 64, 84, 78, 73, 71, 58, 46, 48, 40, 40, 58, 55, 51, 46, 63, 66, 60, 50, 39, 37, 41, 48, 66, 72, 63, 52, 53, 36, 22, 18, 19, 7, 82, 83, 71, 73, 78, 74, 80, 84, 79, 72, 73, 73, 80, 84, 73, 60, 65, 55, 49, 48, 18, 23, 14, 22, 30, 35, 27, 31, 37]]
#divide2 = [[0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]]


finalDivideprint = []
timeX = []
finalDivideprint2 = []
timeX2 = []
for i in range(len(locs[0])):

    if divide[0][i] == 1:
        finalDivideprint.append(locs[0][i])
        timeX.append(i)
    if divide2[0][i] == 1:
        finalDivideprint2.append(locs2[0][i])
        timeX2.append(i)

finalLocs1 = []
finalLocs2 = []
finalx1 = []
finalx2 = []
currlocs1 = []
currlocs2 = []
currx1 = []
currx2 = []
for i in range(len(locs[0])):
    if splits[0][i] == -1:
        currlocs1.append(-50)
        currx1.append(i)
        finalLocs1.append(currlocs1)
        finalx1.append(currx1)
        currx1 = []
        currlocs1 = []
        currlocs1.append(150)
        currlocs1.append(locs[0][i])
        currx1.append(i-1)
        currx1.append(i)

    if splits[0][i] == 1:
        currlocs1.append(150)
        currx1.append(i)
        finalLocs1.append(currlocs1)
        finalx1.append(currx1)
        currx1 = []
        currlocs1 = []
        currlocs1.append(-50)
        currlocs1.append(locs[0][i])
        currx1.append(i - 1)
        currx1.append(i)
    if splits[0][i] == 0:
        currlocs1.append(locs[0][i])
        currx1.append(i)
finalLocs1.append(currlocs1)
finalx1.append(currx1)

for i in range(len(locs2[0])):
    if splits2[0][i] == -1:
        currlocs2.append(-50)
        currx2.append(i)
        finalLocs2.append(currlocs2)
        finalx2.append(currx2)
        currx2 = []
        currlocs2 = []
        currlocs2.append(150)
        currlocs2.append(locs2[0][i])
        currx2.append(i-1)
        currx2.append(i)
    if splits2[0][i] == 1:
        currlocs2.append(150)
        currx2.append(i)
        finalLocs2.append(currlocs2)
        finalx2.append(currx2)
        currx2 = []
        currlocs2 = []
        currlocs2.append(-50)
        currlocs2.append(locs2[0][i])
        currx2.append(i - 1)
        currx2.append(i)
    if splits2[0][i] == 0:
        currlocs2.append(locs2[0][i])
        currx2.append(i)
finalLocs2.append(currlocs2)
finalx2.append(currx2)
print("1-----------------------")
print(finalLocs1)
print()
print(finalx1)
print("2-----------------------")
print(finalLocs1)
print()
print(finalx1)

#plt.scatter(timeX, finalDivideprint)
#plt.scatter(timeX2, finalDivideprint2)
#plt.plot(x, locs[0])
#plt.plot(x, locs2[0])
#plt.show()
#plt.clf()
EnviornmentParams = [100, True, True]
CellMetaStats[0] = 3
simRepeat = 1
allLocations = []
allDivides = []
for i in range(simRepeat):
    print(str(i) + ": Repeat Number")
    timeRun= run.moveTrackRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True, True)
    allLocations = allLocations + timeRun[0]
    allDivides = allDivides + timeRun[1]

print(allLocations)
print(allDivides)

finalDivides = [0]*SimParams[0]
fullLocations = [0] *SimParams[0]
allLocSum = 0
allDivideSum = 0
for i in range(len(allLocations)):
    #print(i)
    #print(allLocations)
    #print(fullLocations)
    fullLocations[allLocations[i]] += 1
    allLocSum += 1
    if allDivides[i] == 1:
        finalDivides[allLocations[i]] = finalDivides[allLocations[i]] + 1
        allDivideSum += 1

finalDivides = [number / allDivideSum for number in finalDivides]
fullLocations = [number / allLocSum for number in fullLocations]

print(finalDivides)

#plt.barh(range(100),finalDivides, height = 1.0, alpha = 0.5)
#plt.barh(range(100),fullLocations, height = 1.0, alpha = 0.5)
#plt.show()

#exit()
fig = plt.figure()
gs = fig.add_gridspec(1, 3, wspace=0,width_ratios=[2, 5, 2])
ax0, ax1, ax2 = gs.subplots( sharey=True)

#fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#fig.add_gridspec(2, wspace=0)
#fig.suptitle('Concentration Distribution Effect on Cell Locations and Division Events (Measured Type)')
print(len(timeRun[11]))
print(len(range(100)))
#concs[0] = concs[0][0:len(concs[0])-1]
#concs[1] = concs[1][0:len(concs[1])-1]
print(len(timeRun[11]))
ax0bar1 = ax0.barh(range(100), timeRun[11], height = 1.0, alpha = 0.5)
ax0bar2 = ax0.barh(range(100), timeRun[12], height = 1.0, alpha = 0.5)
ax0.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax0.set_yticklabels([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], rotation = 270)
ax0.yaxis.set_label_coords(-0.5, 0.5)
ax0.set_ylabel("Location ($i$)", rotation = 270)
ax0.set_xlabel("Concentration", rotation = 270)
ax0.set_xticks([0,5,10,15])
ax0.set_xticklabels([0,5,10,15], rotation = 270)
ax0.xaxis.tick_top()
ax0.xaxis.set_label_position('top')
ax0.legend((ax0bar1, ax0bar2), ('$[A]_i$', '$[B]_i$'), bbox_to_anchor=(7, 0.5))
print()
print(finalDivideprint)

ax1.scatter(timeX, finalDivideprint, marker = 'x', color = 'c')
ax1.scatter(timeX2, finalDivideprint2, marker = 'x', color = 'm')

for i in range(len(finalLocs1)):
    ax1.plot(finalx1[i], finalLocs1[i], color = 'c')

print(finalLocs2)

for i in range(len(finalLocs2)):
    ax1.plot(finalx2[i], finalLocs2[i], color = 'm')

ax1.set_ylim(0, 100)
ax1.grid()

ax1.set_xticks(ticks = [10, 20, 30, 40, 50,  60, 70, 80, 90])
ax1.set_xticklabels([10, 20, 30, 40, 50,  60, 70, 80, 90], rotation = 270)
ax1.xaxis.tick_top()
ax1.set_xlabel("Time Step ($j$)", rotation = 270)
ax1.xaxis.set_label_position('top')
#ax1.legend( ('cell 1 location', 'cell 2 location'), loc='upper right')
#ax1.legend( ('cell 1 locations/divide events', 'cell 2 locations/divide events'), loc='upper left')

ax2plot1 = ax2.barh(range(100),finalDivides, height = 1.0, alpha = 0.5)
ax2plot2 = ax2.barh(range(100),fullLocations, height = 1.0, alpha = 0.5)
ax2.set_xticks(ticks = [0.00, 0.01])
ax2.set_xticklabels([0.00, 0.01], rotation = 270)
ax2.xaxis.tick_top()
ax2.set_xlabel("Cell Event \n Density", rotation = 270)
ax2.xaxis.set_label_position('top')
ax2.legend((ax2plot1, ax2plot2), ('Divide Density', 'Cell Density'), bbox_to_anchor=(1, 0.5))

#plt.show()
plt.rc('pdf', fonttype=42)
fig.savefig('samplefig.pdf',bbox_inches='tight')

