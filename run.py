import simulation
import random
import randomClass
import math
import utils
import itertools
import output
import matplotlib.pyplot as plt
import log
l = log.log()

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

    l.log(Celllocations)

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
    l.log(Celllocations)
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
    l.log(out[1])
    l.log(out[2])

def runIndividual(run, date, config, stress = None, strat = None, envornment = None, Aratio = None, AratioInt = None):
    if stress:
        config.cellMetaStats.absorptionRate = 1 / (stress / config.cellMetaStats.survivalCost)
        config.cellMetaStats.stress = stress
    if strat:
        config.cellMetaStats.decisiontype = strat
    if envornment:
        config.concParams.concProfile = envornment
    if Aratio:
        config.cellMetaStats.Aratio = Aratio
        #config.cellStats.Arec = math.floor(Aratio*config.cellStats.maxRec)
        #config.cellStats.Brec = config.cellStats.maxRec - config.cellStats.Arec
    if AratioInt:
        config.cellMetaStats.AratioInt = AratioInt
    Sim = simulation.simulation(config)
    return(Sim.test(config, run, date))


def testRun(config, runName, date, specificName = None):

    utils.createDirectory(runName, date, config)
    utils.createDirectory(runName + "/OutputProfiles", date, config)
    run_Param_list = []

    if config.runStats.cellStrategiesArrayFlag:
        l.log("Running Strategies: " + str(config.runStats.cellStrategiesArray))
        run_Param_list.append(config.runStats.cellStrategiesArray)
    else:
        run_Param_list.append([None])

    if config.runStats.enviornmentArrayFlag:
        l.log("Running Enviornments: " + str(config.runStats.enviornmentArray))
        run_Param_list.append(config.runStats.enviornmentArray)
    else:
        run_Param_list.append([None])

    if config.runStats.runStress:
        l.log("Running Stress Array: " + str(config.runStats.stressArray))
        run_Param_list.append(config.runStats.stressArray)
    else:
        run_Param_list.append([None])

    if config.runStats.cellRatioAEmphasisFlag:
        l.log("Running A Ratio Array: " + str(config.runStats.cellRatioAEmphasis))
        run_Param_list.append(config.runStats.cellRatioAEmphasis)
    else:
        run_Param_list.append([None])

    if config.runStats.cellRatioAIntEmphasisFlag:
        l.log("Running A Ratio Int Array: " + str(config.runStats.cellRatioAIntEmphasis))
        run_Param_list.append(config.runStats.cellRatioAIntEmphasis)
    else:
        run_Param_list.append([None])

    combinations = itertools.product(*run_Param_list)
    outputFiles = []
    outputObj = []

    count = 0
    for combo in combinations:
        count += 1
        l.log(combo)
        stress = None
        strat = None
        envior = None
        Aratio = None
        AratioInt = None

        print_str = ""
        run_str = ""
        if config.runStats.cellStrategiesArrayFlag:
            strat = combo[0]
            print_str += "Strategy: " + strat + " "
            run_str += strat
        if config.runStats.enviornmentArrayFlag:
            envior = combo[1]
            if len(run_str) > 0:
                run_str += '_'
            print_str += "Enviornment: " + envior + " "
            run_str += envior
        if config.runStats.runStress:
            stress = combo[2]
            if len(run_str) > 0:
                run_str += '_'
            print_str +=  "Stress: " + str(stress) + " "
            run_str += "Stress"+ str(stress)
        if config.runStats.cellRatioAEmphasisFlag:
            Aratio = combo[3]
            if len(run_str) > 0:
                run_str += '_'
            print_str += "Aratio: " + str(Aratio) + " "
            run_str += "Aratio" + str(Aratio)
        if config.runStats.cellRatioAIntEmphasisFlag:
            AratioInt = combo[4]
            if len(run_str) > 0:
                run_str += '_'
            print_str += "AratioInt: " + str(AratioInt) + " "
            run_str += "AratioInt" + str(AratioInt)
        l.log(print_str)
        run = runName + "/" + run_str + "_" + utils.getTime()
        if specificName:
            run += "_" + specificName
        runResults = runIndividual(run, date, config, stress=stress, strat=strat, envornment=envior, Aratio=Aratio, AratioInt=AratioInt)
        outputObj.append(runResults[0])
        outputFiles.append(runResults[1])
    return outputObj, outputFiles

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

    l.log(Celllocations)

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

    l.log(Celllocations)

    MIs = []
    l.log(CellStats[2])
    for i in range(5, CellStats[2],10):
        l.log("A receptors: " + str(i))
        CellStats[0] = i
        CellStats[1] = CellStats[2] - CellStats[0]
        l.log(CellStats[0])
        l.log(CellStats[1])
        Sim = simulation.simulation(ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams,
                                SimParams, newRand)
        div, MI = Sim.staicConcRunTime(flag)
        MIs.append(MI)
    return [MIs, randomSeed]