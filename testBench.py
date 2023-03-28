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
import random
import log

l = log.log()

def main():

    meta = utils.loadMetaData()
    givenRunName = meta["givenRunName"]
    runSim = meta["runSim"]
    newRun = meta["newRun"]
    refactorPaths = meta["refactorOutputProfilePaths"]
    analyzeResults = meta["analyzeResults"]
    givenDate = meta["givenDate"]
    calculateResults = meta["calculateResults"]
    calculation_parameters = meta["calculation_parameters"]
    config_file = meta["config_file"]

    l.log("Run Name: " + givenRunName + " Date: " + givenDate)
    if runSim:
        l.log('Running ' + '\'' + givenRunName+ '\'')
    if analyzeResults:
        l.log('Analyzing '+ '\'' + givenRunName+ '\'')
    if calculateResults:
        l.log('Calculating Results '+ '\'' + givenRunName+ '\'')

    c = configuration.configuration()
    c.readConfig(config_file)
    c.runStats.saveDir = meta['path']

    # give run name
    runName = givenRunName
    date = givenDate

    output_objects = []
    output_files = []
    if runSim:
        if newRun:
            date = str(utils.getTodaysDate())
        output_objects, output_files = run.testRun(c, runName, date)
        #utils.saveData(c.runStats.saveDir +"/" + date + "/" + runName + "/OutputProfiles" , output_files, "outputProfilePaths.txt")

    if not utils.pathExists(c.runStats.saveDir + '/' + date + '/' + runName):
        l.log("Given run Path: " + "\'"+c.runStats.saveDir + '/' + date + '/' + runName+ "\' does not exist. May need to run the simulation")
        runNames = utils.findRunNames(c.runStats.saveDir)
        if not runNames == []:
            l.log("Here are some Runs found on the System: ")
            l.log(runNames)
            return

    output_files = utils.loadPlain(c.runStats.saveDir + "/" + date + "/" + runName + "/OutputProfiles/outputProfilePaths.txt")
    for i in range(len(output_files)):
        file = (c.runStats.saveDir + '/' + output_files[i]).strip()
        output_files[i] = file
        o = output.output()
        o.read(file)
        output_objects.append(o)

    if calculateResults:
        #if not newRun or not runSim:
        output_files = utils.loadPlain(c.runStats.saveDir + "/" + date + "/" + runName + "/OutputProfiles/outputProfilePaths.txt")
        for i in range(len(output_files)):
            file = (c.runStats.saveDir + '/' + output_files[i]).strip()
            output_files[i] = file
            o = output.output()
            o.read(file)
            output_objects.append(o)
        l.log("Calculating Results on Parameters run")
        calcParamDict = {0 : 'Traditional MI Calculation', 1 : 'New Syntactic Caclulation', 2 : 'New Useful Calculation', 3: 'New Useful Input Calculation', 4 :  'Growth Calculation', 5 : 'New Internal Weighted Calculation', 6: 'Dynamic Weight'}
        ret = 'Starting: '
        for i in range(len(calculation_parameters)):
            if calculation_parameters[i]:
                ret += calcParamDict[i] + ', '
        l.log(ret)
        bins = 30
        for i, out in enumerate(output_objects):
            l.log("calculating...")
            out.calculateMeasures(c, bins, MItradFlag = calculation_parameters[0], MI2D2DFlag = calculation_parameters[1], MI2D1DFlag = calculation_parameters[2], MI2D1DinFlag = calculation_parameters[3], growthFlag = calculation_parameters[4], intweightMImoveFlag = calculation_parameters[5], ABintdynamicWeightFlag = calculation_parameters[6], path = c.runStats.saveDir)
            l.log(out.write(output_files[i], absolute=True))

    if analyzeResults:
        outputAnalysis.createFigureGrowthMISyntactic(c, output_objects)
    """""

    #plt.show()

    l.exit()
    plt.scatter(MImoves, growths, color=colors_plot2)
    plt.scatter(2, 2, color='Blue', label="Internal Ratio")
    plt.scatter(2, 2, color='Red', label="Receptor Ratio")
    plt.xlabel("MI (2D1D)")
    # plt.xlim([0, 0.6])
    plt.ylabel("growth")
    # plt.ylim([-.1, .25])
    plt.grid()
    plt.legend()
    plt.show()

    return
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
    """""
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

def createNewConfig(configName):
    c_presets = configuration.configuration()

    example_configuration = "Data/2023-02-06/manaFromHeavenPresets/single/config.cfg"
    # load in example config_file
    c_presets.readConfig(example_configuration)

    # create new default configuration file
    c = configuration.configuration()

    # set mana From Heaven Presets (see README.txt for more information)
    c.simParams.presetA_time = c_presets.simParams.presetA_time
    c.simParams.presetB_time = c_presets.simParams.presetB_time
    c.simParams.presetA_loc = c_presets.simParams.presetA_loc
    c.simParams.presetB_loc = c_presets.simParams.presetB_loc

    # set config flag to save results
    c.runOutputFlags.save = True
    # set simulation time length
    c.simParams.simLength = 10
    # set save compression
    c.runOutputFlags.compressSave = True
    # set simulation time step
    c.simParams.simTimeStep = .05
    c.cellMetaStats.stress = 0.05
    c.cellMetaStats.absorptionRate = 1 / c.cellMetaStats.stress
    c.cellMetaStats.survivalCost = 100
    c.concParams.VonMisesMagnitude = 200
    c.runOutputFlags.compressSave = False
    # set Cell strategy
    c.cellMetaStats.decisiontype = "non"

    c.simParams.highSpace = 200
    c.simParams.lowSpace = 100

    # set meta simulation parameters
    c.runStats.runStress = False
    c.runStats.cellStrategiesArrayFlag = True
    c.runStats.enviornmentArrayFlag = False
    c.runStats.cellRatioAEmphasisFlag = True
    c.runStats.cellRatioAIntEmphasisFlag = True
    # TODO: implement below
    c.runStats.runDetermine = False

    # set other simulation parameters
    c.simParams.presetSave = False
    c.concParams.concProfile = 'vonMises'
    c.simParams.presetBool = True
    c.cellMetaStats.dissocociationConstant = 2.0

    # set output file flags
    c.outputFlags.totalcells = True
    c.outputFlags.totalConcs = True
    c.outputFlags.totalReceptors = True
    c.outputFlags.boundReceptors = True
    c.outputFlags.internalAB = True
    c.outputFlags.cellLocations = True
    c.outputFlags.cellMovement = True

    # set simulation meta parameters
    c.runStats.stressArray = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #c.runStats.cellRatioAEmphasis = [-50, -40, -30, -20, -10, -5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1,1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 20, 30, 40, 50]
    c.runStats.cellRatioAEmphasis = [-50,-40, -30, -20, -10, -5, -2, -1, 0, 1,2,5, 10, 20, 30, 40, 50]
    #c.runStats.cellRatioAEmphasis = [-11,-9, -7]
    #c.runStats.cellRatioAIntEmphasis = [-50, -40, -30, -20, -10, -5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1,1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 20, 30, 40, 50]
    c.runStats.cellRatioAIntEmphasis = [-50,-40, -30, -20, -10, -5, -2, -1, 0, 1,2,5, 10, 20, 30, 40, 50]
    #c.runStats.cellRatioAIntEmphasis = [-11,-9, -7]
    c.runStats.enviornmentArray = ["manaFromHeaven"]
    c.runStats.cellStrategiesArray = ["ratio"]
    path = 'Configs/'+configName
    c.writeConfig(path)
    return path