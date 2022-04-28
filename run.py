import simulation
import random
import randomClass

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

    div, MI, rolling, entr_in, hx, hxgiveny, MI_vars = Sim.staicConcRunTimeVar(True)
    return [div, MI, rolling, entr_in,hx, hxgiveny,  randomSeed, MI_vars]

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