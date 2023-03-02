"""""
#testBench.py
#Tyler Barker
#2-10-2023
"""""
import typing

import utils
import configuration
import run
import outputAnalysis
import output
import math
from scipy import special
import matplotlib.pyplot as plt
import numpy as np
import colors
import log

dateTime = str(utils.getTodaysDate()) + '_' + utils.getTime()
l = log.log('Logs/'+dateTime + '.log')

"""""
#----------------------------
finalDat = utils.loadDataDate("Data/finalDataSave_08_04_00.834293.gz", True)[0]
#growths = [[ 1.3, 1.2, 1.1, 1.4, 1.0, .7, .65, .55, .42, .35], [ 1.4, 1.3, 1.2, 1.1, 1.0, .9, .7, .6, .5, .2], [ .9, .95, 1.0, .87, .80, .7, .5, .35, .3, .27], [ .6, .65, .73, .80, .9, .75, .6, .63, .5, .45]]
#MIs = [[ 1, 1.1, 1.4, 1.6, 1.8, 1.9, 1.5, 1.7, 1.9, 3], [ .6, .89, .9, .7, .85, .95, 1.5, 1.75, 2.5, 4], [ .6, .65, .73, .80, .9, .75, .6, .63, .5, .45], [ .9, .95, 1.0, .87, .80, .7, .5, .35, .3, .27]]
stresses = [0.05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
MI_trad = np.array(finalDat[0])[:,2,:]
print(len(MI_trad))
MImove = np.array(finalDat[2])[:,2,:]
growths = np.array(finalDat[4])[:,2,:]
strategies = ["Equal", "Adaptive", "Adjusted", "Drastic", "Drastic2"]
#outputAnalysis.createFigureMIGrowth(growths, MImove, MImove, strategies, stresses)
for str in stresses:
    outputAnalysis.createFigureMIGrowthStress(growths, MImove, MImove, strategies, stresses, str)
exit()
#----------------------------
"""""

"""""
def VonMises(var, maximum, offset, length, locationStep):
    arr = []
    kappa = 1 / var
    mu = (2 * math.pi / length) * offset
    bessel = (special.iv(0, kappa))
    sum = 0.0
    increment = 2 * math.pi / length
    for i in range(length):
        index = i - length / 2
        x = index * increment
        numerator = math.exp(kappa * math.cos(x - mu))
        arr.append((numerator / (2 * math.pi * bessel)) * maximum)
        sum += numerator / (2 * math.pi * bessel)
    #arr = (arr / sum) * maximum
    return arr

var = 10 * ((2*math.pi)/100)
maximum = 10
#maximum = 10
locationStep = 0.1
offset1 = -25 * (1/locationStep)
offset2 = 25 * (1/locationStep)
length = int(100 * (1/locationStep))
arr1 = VonMises(var, maximum, offset1, length, locationStep)
arr2 = VonMises(var, maximum, offset2, length, locationStep)

plt.plot(arr1)
plt.plot(arr2)

minimumi = 0
minimum = 100

for i in range(len(arr1)):
    if math.fabs(arr2[i] - arr1[i]) <  minimum:
        minimumi = i
        minimum = math.fabs(arr2[i] - arr1[i])

print(arr1[minimumi])
equal = arr1[minimumi]

stress_rate = 0.9
timeDivison = 0.01
survival_cost = 250
absoption_rate = 1/stress_rate
adjustedcost = timeDivison*survival_cost
adjusted_absorption_rate = absoption_rate * timeDivison
print("adjusted_absorption_rate: " + str(adjusted_absorption_rate))
absorbed = (adjusted_absorption_rate * equal)*(1/timeDivison)
print("cost: " + str(adjustedcost))
print("absorbed: " + str(absorbed))

plt.show()



exit()
maximum = 10
stress_rate = 0.9
"""""

c_presets = configuration.configuration()

example_configuration = "Data/2023-02-06/manaFromHeavenPresets/single/config.cfg"
#load in example config_file
c_presets.readConfig(example_configuration)

#create new default configuration file
c = configuration.configuration()

#set mana From Heaven Presets (see README.txt for more information)
c.simParams.presetA_time = c_presets.simParams.presetA_time
c.simParams.presetB_time = c_presets.simParams.presetB_time
c.simParams.presetA_loc = c_presets.simParams.presetA_loc
c.simParams.presetB_loc = c_presets.simParams.presetB_loc

#set config flag to save results
c.runOutputFlags.save = True
#set simulation time length
c.simParams.simLength = 20
#set save compression
c.runOutputFlags.compressSave = True
# set simulation time step
c.simParams.simTimeStep = .05
c.cellMetaStats.stress = 0.05
c.cellMetaStats.absorptionRate = 1/c.cellMetaStats.stress
c.cellMetaStats.survivalCost = 100
c.concParams.VonMisesMagnitude = 200
c.runOutputFlags.compressSave = False
# set Cell strategy
c.cellMetaStats.decisiontype = "non"

# set meta simulation parameters
c.runStats.runStress = False
c.runStats.cellStrategiesArrayFlag = True
c.runStats.enviornmentArrayFlag = False
c.runStats.cellRatioAEmphasisFlag = False
c.runStats.cellRatioAIntEmphasisFlag = False
# TODO: implement below
c.runStats.runDetermine = False

#set other simulation parameters
c.simParams.presetSave = False
c.concParams.concProfile = 'vonMises'
c.simParams.presetBool = True
c.cellMetaStats.dissocociationConstant = 2.0

#set output file flags
c.outputFlags.totalcells = True
c.outputFlags.totalConcs = True
c.outputFlags.totalReceptors = True
c.outputFlags.boundReceptors = True
c.outputFlags.internalAB = True
c.outputFlags.cellLocations = True
c.outputFlags.cellMovement = True

#set simulation meta parameters
c.runStats.stressArray = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
c.runStats.enviornmentArray = ["vonMises"]
c.runStats.cellStrategiesArray = ["ratio"]


meta = utils.loadMetaData()

c.runStats.saveDir = meta["path"]

# give run name
runName = 'output_save_test'
date = str(utils.getTodaysDate())

#output_objects, output_files = run.testRun(c, runName, date)
#utils.saveData(c.runStats.saveDir +"/" + date + "/" + runName + "/OutputProfiles" , output_files, "outputProfilePaths.txt")

date = "2023-02-28"
runName = "ratioA_IntA"
output_files = utils.loadData(c.runStats.saveDir +"/"+ date+ "/" + runName  + "/OutputProfiles/outputProfilePaths.txt")


#TODO needs to be rebuilt
def refactorfile(file):
    arr = file.split('/')
    count = 0
    for a in arr:
        count += 1
        if a == 'Data':
            break
    newFile = ''
    for i in range(count, len(arr), 1):
        newFile += arr[i]
        if i != len(arr)-1:
            newFile += '/'
    return newFile

output_objects = []
for i in range(len(output_files)):
    file = c.runStats.saveDir + '/' + refactorfile(output_files[i])
    l.log(file)
    output_files[i] = file
    o = output.output()
    o.read(file)
    output_objects.append(o)

"""""
bins = 30
for i, out in enumerate(output_objects):
    out.calculateMeasures(c, bins)
    print("finished: " + out.runName)
    print(out.write(output_files[i], absolute=True))
    #exit()
    
"""""
def d3color(emp_red, emp_blue):
    color_Scheme1 = ['Black', 'Red']
    color_Scheme2 = ['Black', 'Blue']
    color_red = colors.findcolor(1,0,color_Scheme1, emp_red, absolute=False)
    color_blue = colors.findcolor(1,0,color_Scheme2, emp_blue, absolute=False)
    color_Scheme3 = [color_red, color_blue]
    dist = color_blue[2] / (color_red[0] + color_blue[2])
    color_purp = colors.findcolor(1,0,color_Scheme3, dist)
    return color_purp


#ratioA is red
#ratioint is blue
colors_plot = []
MImoves = []
MIs = []
MITrad = []
growths = []
count = -1
for i in range(len(c.runStats.cellRatioAEmphasis)):
    for j in range(len(c.runStats.cellRatioAIntEmphasis)):
        count += 1
        if output_objects[count].RunClacOut.MI2D2D != '':
            colors_plot.append(d3color(c.runStats.cellRatioAEmphasis[i], c.runStats.cellRatioAIntEmphasis[j]))
            MImoves.append(output_objects[count].RunClacOut.MI2D1D)
            MIs.append(output_objects[count].RunClacOut.MI2D2D)
            MITrad.append(output_objects[count].RunClacOut.MItrad)
            growths.append(output_objects[count].RunClacOut.growth)

plt.scatter(MIs, growths, color = colors_plot)
plt.scatter(2, 2, color = 'Blue', label = "Internal Ratio")
plt.scatter(2, 2, color = 'Red', label = "Receptor Ratio")
plt.xlabel("MI (2D2D)")
plt.xlim([0, 0.6])
plt.ylabel("growth")
plt.ylim([-.1, .25])
plt.grid()
plt.legend()
plt.show()


l.exit()
bins = 30
stresses = c.runStats.stressArray
strats = c.runStats.cellStrategiesArray
enviorns = c.runStats.enviornmentArray

#analyses and saves data (see README.txt for more imformation)
finalDat = outputAnalysis.getCalcArray(c, runName, date, stresses, strats, enviorns, bins)
utils.saveDataDate(runName, date,"finalDataSave", finalDat, True, c.runStats.saveDir)

# can load in particular data
#finalDat = utils.loadDataDate("Data/finalDataSave_08_04_00.834293.gz", True)[0]

#creates figures for data in \Results
outputAnalysis.createFigures(strats, stresses, enviorns, finalDat, runName)

# uncomment below to create gif of enviornment and cells over time (example gif in \sim_gifs)
"""""
enviornmentFile = "Data/2023-02-06/manaFromHeavenTest2/single/enviornmentConcs_17_25_24.042995.gz"
cellsFile = "Data/2023-02-06/manaFromHeavenTest2/single/cellLocations_17_25_23.527882.gz"

enviornmentConcs = utils.loadDataDate(enviornmentFile, True)
cellslocs = utils.loadDataDate(cellsFile, True)[0]

length = int(c.simParams.length * (1/c.simParams.locationStep))

length_time = len(cellslocs)
length_space = length
cellArr = utils.list_empty([length_time, length_space], 0)
for i in range(len(cellslocs)):
    for j in range(len(cellslocs[i])):
        cellArr[i][cellslocs[i][j]] += 1

utils.createGifBar(c, enviornmentConcs[0], cellArr, runName)
"""""