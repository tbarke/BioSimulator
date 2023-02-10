import simulation
import random
import randomClass
import math
import utils
import output
import matplotlib.pyplot as plt

def moveTrackRunReceptorIntAB(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag, flag2):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    Celllocations = []
    if not flag2 :
        Celllocations = [51]
    if flag2:
        for j in range(10):
            for i in range(SimParams[0]):
                Celllocations.append(i)
    CellLocationStats = [Celllocations]

    print(Celllocations)

    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
    ret = Sim.moveRunReceptors()
    arr = [Sim.SimEnviornment.Aconcentrations, Sim.SimEnviornment.Bconcentrations, Sim.SimEnviornment.IntA, Sim.SimEnviornment.IntB, Sim.SimEnviornment.divideLoc]
    return [ret , arr]

def moveTrackRunReceptor(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag, flag2):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    Celllocations = []
    if not flag2 :
        Celllocations = [51]
    if flag2:
        for j in range(10):
            for i in range(SimParams[0]):
                Celllocations.append(i)
    print(Celllocations)
    exit(-1)
    CellLocationStats = [Celllocations]

    #print(Celllocations)

    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
    concs = [Sim.SimEnviornment.Aconcentrations, Sim.SimEnviornment.Bconcentrations]
    return [Sim.moveRunReceptors(), concs]

def moveTrackRunAB(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag, flag2):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    Celllocations = []
    if not flag2 :
        Celllocations = [51]
    if flag2:
        for j in range(1):
            for i in range(SimParams[0]):
                Celllocations.append(i)
    CellLocationStats = [Celllocations]

    #print(Celllocations)

    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
    return Sim.moveRunAB()

def moveTrackRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag, flag2):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    Celllocations = []
    if not flag2 :
        Celllocations = [51]
    if flag2:
        for j in range(1):
            for i in range(SimParams[0]):
                Celllocations.append(i)
    CellLocationStats = [Celllocations]

    print(Celllocations)

    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
    concs = [Sim.SimEnviornment.Aconcentrations, Sim.SimEnviornment.Bconcentrations]
    return Sim.moveRun()

def indRun(config):

    run = "k"

    #intB = 500-700
    #intA = 400 - 600
    config.outputFlags.boundReceptors = False
    config.outputFlags.cellMovement = False
    config.outputFlags.totalConcs = False
    config.outputFlags.totalcells = False
    config.outputFlags.totalReceptors = False
    config.outputFlags.cellLocations = False
    stress = .1
    config.cellMetaStats.absorptionRate = 1 / (stress / config.cellMetaStats.survivalCost)
    config.cellMetaStats.stress = stress
    config.cellMetaStats.decisiontype = "measured"
    Sim = simulation.simulation(config)
    out = Sim.test(config, run)
    intABFile = out.primativeOutput.internalAB
    intAB = utils.loadDataDate(intABFile, config.runOutputFlags.compressSave, not config.runOutputFlags.save,
                             not config.runOutputFlags.save)[0]

    locs_A = []
    #print(intAB)
    for i in range(len(intAB)):
        for j in range(len(intAB[i])):

            if intAB[i][j][0] > 500 and intAB[i][j][0] < 700 and intAB[i][j][1] > 400 and intAB[i][j][1] < 600:
                locs_A.append(j)

    config.cellMetaStats.cellLocations = locs_A

    #plt.hist(locs_A, bins = 50)
    #plt.show()
    #exit(-1)

    config.cellMetaStats.fullDivide = False
    Sim = simulation.simulation(config)
    out = Sim.indtest(config, run)
    print(out[1])
    print(out[2])

def runIndividual(run, date, config, stress = None, strat = None, envornment = None):
    if stress:
        config.cellMetaStats.absorptionRate = 1 / (stress / config.cellMetaStats.survivalCost)
        config.cellMetaStats.stress = stress
    if strat:
        config.cellMetaStats.decisiontype = strat
    if envornment:
        config.concParams.concProfile = envornment
    Sim = simulation.simulation(config)
    return(Sim.test(config, run, date))


def testRun(config, runName, date):
    #stress_arr = [.01, .02, .03, .04, .05, .06, .07, .08, .09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2.0]
    #stress_arr = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    #cell_strategies = ["non", "measured", "counter", "adjusted", "drastic"]
    #cell_strategies = ["drastic2"]
    #TODO
    #config.runStats.runNoise

    date = utils.createDirectory(runName)

    if config.runStats.cellStrategiesArrayFlag:
        print("Running Strategies: " + str(config.runStats.cellStrategiesArray))
    if config.runStats.enviornmentArrayFlag:
        print("Running Enviornments: " + str(config.runStats.enviornmentArray))
    if config.runStats.runStress:
        print("Running Stress Array: " + str(config.runStats.stressArray))

    #need a better way to do this
    ouputFileObjects = []
    if config.runStats.cellStrategiesArrayFlag:
        if config.runStats.enviornmentArrayFlag:
            if config.runStats.runStress:
                #running all three
                for strat in config.runStats.cellStrategiesArray:
                    for envior in config.runStats.enviornmentArray:
                        for stress in config.runStats.stressArray:
                            print("Stress: " + str(stress) + " Strategy: " + strat + " Enviornment: " + envior)
                            run = runName + "/" + str(stress) + "_" + strat + "_" + envior
                            ouputFileObjects.append(runIndividual(run, date, config, stress=stress, strat=strat, envornment= envior))
            else:
                #run strategy and enviornment
                for strat in config.runStats.cellStrategiesArray:
                    for envior in config.runStats.enviornmentArray:
                        print("Strategy: " + strat + " Enviornment: " + envior)
                        run = runName + "/" + strat + "_" + envior
                        ouputFileObjects.append(runIndividual(run, date, config, strat=strat, envornment=envior))
        elif config.runStats.runStress:
            for strat in config.runStats.cellStrategiesArray:
                for stress in config.runStats.stressArray:
                    print("Stress: " + str(stress) + " Strategy: " + strat)
                    run = runName + "/" + str(stress) + "_" + strat
                    ouputFileObjects.append(runIndividual(run, date, config, stress=stress, strat=strat))
        else:
            #run strategy
            for strat in config.runStats.cellStrategiesArray:
                print("Strategy: " + strat)
                run = runName + "/" + strat
                ouputFileObjects.append(runIndividual(run, date, config, strat=strat))

    elif config.runStats.enviornmentArrayFlag:
        if config.runStats.runStress:
            #run strategy enviornment and stress
            for envior in config.runStats.enviornmentArray:
                for stress in config.runStats.stressArray:
                    print("Stress: " + str(stress) + " Enviornment: " + envior)
                    run = runName + "/" + str(stress)+ "_" + envior
                    ouputFileObjects.append(runIndividual(run, date, config, stress=stress, envornment=envior))
        else:
            #run enviornment
            for envior in config.runStats.enviornmentArray:
                print("Stress: "  + " Enviornment: " + envior)
                run = runName + "/" + envior
                ouputFileObjects.append(runIndividual(run, date, config, envornment=envior))
    elif config.runStats.runStress:
        #run stress
        for stress in config.runStats.stressArray:
            print("Stress: " + str(stress))
            run = runName + "/" + str(stress)
            ouputFileObjects.append(runIndividual(run, date, config, stress=stress))
    else:
        #single run
        Sim = simulation.simulation(config)
        ouputFileObjects.append(Sim.test(config, runName + "/" + "single", date))

    return ouputFileObjects

def testRunConcs(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    Celllocations = []
    for j in range(1):
        for i in range(math.floor(SimParams[0]/EnviornmentParams[4])):
            Celllocations.append(i)
    #print(Celllocations)

    # Celllocations.append(50)
    CellLocationStats = [Celllocations]
    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)

    return Sim.testconc()

def divisonMIRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    #Celllocations = [51]
    Celllocations = []
    for j in range(1):
        for i in range(SimParams[0]):
            Celllocations.append(i)

    #Celllocations.append(50)
    CellLocationStats = [Celllocations]

    #print(Celllocations)

    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)

    div, MI, rolling, entr_in, hx, hxgiveny, MI_vars, enviornment_Stats = Sim.staicConcRunTimeVar(True)
    return [div, MI, rolling, entr_in,hx, hxgiveny,  randomSeed, MI_vars, enviornment_Stats]

def SNRrun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag):
    randomSeed = random.randint(1, 10000)
    newRand = randomClass.randomClass(randomSeed)
    for j in range(1):
        for i in range(SimParams[0]):
            Celllocations.append(i)
    CellLocationStats = [Celllocations]

    print(Celllocations)

    Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
    div, MI, SNR = Sim.staicConcRunTimeSNR(flag)
    return [div, MI, SNR, randomSeed]

def MI_sensitivity_curve(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, flag):
    randomSeed = random.randint(1, 10000)

    newRand = randomClass.randomClass(randomSeed)
    for j in range(30):
        for i in range(SimParams[0]):
            Celllocations.append(i)
    CellLocationStats = [Celllocations]

    print(Celllocations)

    MIs = []
    print(CellStats[2])
    for i in range(5, CellStats[2],10):
        print("A receptors: " + str(i))
        CellStats[0] = i
        CellStats[1] = CellStats[2] - CellStats[0]
        print(CellStats[0])
        print(CellStats[1])
        Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
        div, MI = Sim.staicConcRunTime(flag)
        MIs.append(MI)
    return [MIs, randomSeed]