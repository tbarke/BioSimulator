"""""
#testBench.py
#Tyler Barker
#2-10-2023
"""""

import utils
import configuration
import run
import outputAnalysis

c_presets = configuration.configuration()
t = " h "
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
# set Cell strategy
c.cellMetaStats.decisiontype = "measured"

# set meta simulation parameters
c.runStats.runStress = True
c.runStats.cellStrategiesArrayFlag = True
c.runStats.enviornmentArrayFlag = True
# TODO: implement below
c.runStats.runDetermine = False

#set other simulation parameters
c.simParams.presetSave = False
c.concParams.concProfile = 'manaFromHeaven'
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
c.runStats.cellStrategiesArray = ["non", "measured", "adjusted", "drastic", "drastic2"]

# give run name
runName = 'vonMises_dis5'
date = str(utils.getTodaysDate())

#outputs output list type containing output files
output_files = run.testRun(c, runName, date)
bins = 30
stresses = c.runStats.stressArray
strats = c.runStats.cellStrategiesArray
enviorns = c.runStats.enviornmentArray

#analyses and saves data (see README.txt for more imformation)
finalDat = outputAnalysis.getCalcArray(c, runName, date, stresses, strats, enviorns, bins)
utils.saveDataDate(runName, date,"finalDataSave", finalDat, True)

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