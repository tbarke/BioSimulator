"""""
#Run the main Simulation here
#Tyler Barker
#3-24-2021
#Demo.py

NEEDED DIRECTORIES under src:
Cells
Concentrations
Output
Output/Data
Output/Plots
Receptors
Receptors/A
Receptors/B
"""""

#from scipy import special
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
import pickle
#import numpy as np
import gc
import run
import time
import runSuite
plt.rc('pdf', fonttype=42)
locationStep = 1
# size scaled from 100
scale = 1/locationStep
SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA, locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25.5 * scale, 25.5 * scale, 0, 50]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change), noise(as a protion of the binomial distribution), receptor_mode(simulate, gaussian), scale, combined_portion, ABdivide_equal]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed)


cell_stress = 1

CellMetaStats = [1.0/cell_stress, 0.0, 1, 1.0, "non", "non", 0, "gaussian", 1, 1, 1]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die, sim Time Step, sim location step]
EnviornmentParams = [50, True, True, 1, locationStep]

CellMetaStats[5] = "non"
simLength = EnviornmentParams[0]
simRepeat = 1
parameters = [SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, simRepeat]

noise_array = [0.0]
cell_stress_array = [0.2, 0.4, 0.6, 0.8, 1]


path = runSuite.suite_noise_stress(noise_array, cell_stress_array, parameters, "Time_Step_Test_final")
data,path = utils.loadDataDate(path, True)
print(path)