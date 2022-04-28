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
from matplotlib.colors import LogNorm
import math
import numpy as np
import gc
import run

scale = 1
SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA (negative if random), locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25.5, 25.5, 0, 50]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change)]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]

CellMetaStats = [5, 0.0, 5, 1.0, "non", "non", 0, "simulate", scale]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []
CellLocationStats = [Celllocations]

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die]
EnviornmentParams = [30, True, True]

CellMetaStats[5] = "measured"
simLength = EnviornmentParams[0]
simRepeat = 1
fianlTimeDivx = range(simLength)
finalTimeDiv = [0]*simLength
survivaly = []

Amol = []
Bmol = []

locs = []
divides = []
splitLoc = []
concs = []
total_Cells = []
divides_per_time = []
mol_times = []
real_cells = []
Acons = []
Bcons = []
for i in range(simRepeat):
    print(str(i) + ": Repeat Number")
    init = run.moveTrackRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, True, True)
    locs = init[0]
    divides = init[1]
    splitLoc = init[2]
    Amol = init[3]
    Bmol = init[4]
    Acons = init[5]
    Bcons = init[6]
    total_Cells = init[7]
    divides_per_time = init[8]
    mol_times = init[9]
    real_cells = init[10]

print(len(Amol))
print(len(Bmol))
print(len(mol_times))
print(len(total_Cells))
print(len(real_cells))
#print(divides)

#print(total_Cells)
#print(divides_per_time)

cell_totals = []
A_totals = []
B_totals = []
for i in range(len(Amol)):
    cell_totals.append(total_Cells[mol_times[i]]/real_cells[mol_times[i]])
    A_totals.append(total_Cells[mol_times[i]]*Acons[i])
    B_totals.append(total_Cells[mol_times[i]]*Bcons[i])

print("amol")
print(Amol)
print("bmol")
print(Bmol)
print("cell_totals")
print(cell_totals)
print("A_totals")

print(A_totals)
print("B_totals")
print(B_totals)
plt.rc('pdf', fonttype=42)
color_scheme = ['p1', 'g1', 'y1']
max = utils.heatMap(Amol, Bmol, cell_totals, 50, 100, 0.004077675082354259, 0, True, "Internal A Molecule Count", "Internal B molecule count", "Cell Density",color_scheme)
print(max)
color_scheme = ['p1', 'g1', 'y1']
max = utils.heatMap(Amol, Bmol, A_totals, 50, 100, 0.003653546496093295, 0, True, "Internal A Molecule Count", "Internal B molecule count", "External A Concentration Density",color_scheme)
color_scheme = ['p1', 'g1', 'y1']
max = utils.heatMap(Amol, Bmol, B_totals, 50, 100, 0.003653546496093295, 0, True, "Internal A Molecule Count", "Internal B Molecule Count", "External B Concentration Density",color_scheme)
print(max)

exit(-1)
newAmol = []
newBmol = []
new_cell_totals = []
total = 0.0
for i in range(len(cell_totals)):
    if Amol[i] < 100 and Bmol[i] < 100:
        newAmol.append(Amol[i])
        newBmol.append(Bmol[i])
        new_cell_totals.append(cell_totals[i])
    total += cell_totals[i]

for i in range(len(new_cell_totals)):
    new_cell_totals[i] = new_cell_totals[i]/total

histAmol = []
histBmol = []

for i in range(len(newAmol)):
    for j in range(int(new_cell_totals[i]*10000)):
        histAmol.append(newAmol[i])
        histBmol.append(newBmol[i])

fig, ax = plt.subplots()
h = plt.hist2d(histAmol,histBmol,bins = 25)
cbar = fig.colorbar(h[3], ax=ax)
plt.xlim(0,100)
plt.ylim(0,100)
plt.show()