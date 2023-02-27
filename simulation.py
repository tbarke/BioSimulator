import cell
import enviornment
import MICalc
import utils
import random
import math
from scipy import special
import output
import matplotlib.pyplot as plt

class simulation(object):
    #concnetrations parameter array
    # [diffCoeff, absorb/reflect, repeatFrequency, magnitude, locationA (negative if random), locationB]
        #concentrations diffusion coefficent (0 if no diffuse)
        #conentration boundary conditions (absorb, reflect)
        #concentration degradation or depletion rate due to cells consumption

    #cells parameter array
    #cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, ID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]

        #cell initilization stats ( A/B rec/mol, ATP, biomass threshold)
        #cell locations and counts
        #cell MetaParameters (receptor binding rate, absorbtion rate, consumption rate, costs...)

    #Enviornemnt paramter array: [simulationLengthTime]
        #simulation length

    #simulation parameters: length, data to record (MI data: input concentrations, resulting velocities) etc..

    #...

    def VonMises(self, var, maximum, offset, length, locationStep):
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
        return arr

    #def uniform(self, ):

    #reads in and returns an SNR value for a given set of receptors and concentrations
    def findSNRfromTable(self, receptors, concentrations):
        SNRs = []
        g = open("concentrationsSNR.txt", "r")
        allconcentrations = g.readlines()[0].split(",")
        for i in range(len(allconcentrations)):
            allconcentrations[i] = float(allconcentrations[i])
        g.close()
        f = open("allSNRtable.txt", "r")
        lines = f.readlines()
        f.close()
        print(receptors)
        #print(concentrations)
        for i in range(len(receptors)):
            if i == 1:
                print(self.SimEnviornment.cellsAlive())
                print("-------------------------------------------------------------------------------")
            if receptors[i] == 0:
                SNRs.append(0.0)
                continue
            for j in range(0, len(lines)):
                arrline = lines[j].split(",")
                if int(arrline[0]) == int(receptors[i]):
                    for k in range(len(allconcentrations)):
                        if math.floor(10000*allconcentrations[k]) == math.floor(10000*concentrations[i]):
                            SNRs.append(float(arrline[k+1]))
                            break
                    break
        print(SNRs)
        sum = 0.0
        for i in SNRs:
            sum += i
        return (sum/len(SNRs))

    #initalize the simulation
    def constantConc(self, conc, length):
        arr = []
        for i in range(length):
            arr.append(conc)
        return arr

    def constantUpDown(self, length, low, high, offset):
        rate = (high-low)/(length/2)
        arr = []
        firstUp = int(length/2)
        down = length - firstUp
        for i in range(int(length/2)):
            arr.append(low+(rate*i))
        for i in range(down):
            arr.append(high-(rate*i))
        arr2 = utils.list_empty([len(arr)], 0)
        offsetAdjusted = math.floor((offset/100)*length)
        for i in range(length):
            arr2[i] = arr[(i+offsetAdjusted) % length]
        return arr2

    def __init__(self, config):
        self.config = config
        self.length = int(config.simParams.length * (1/config.simParams.locationStep))
        absorb = None
        if config.concParams.arp == "absorb":
            absorb = True
        self.concentrationFreq = config.concParams.repeatFrequency
        self.concentrationMag = config.concParams.magnitude
        self.locationA = config.concParams.locationA
        self.locationB = config.concParams.locationB
        self.locationC = config.concParams.locationC
        ConcentrationA = []
        ConcentrationB = []
        ConcentrationC = []
        for i in range(self.length):
            ConcentrationA.append(0.0)
            ConcentrationB.append(0.0)
            ConcentrationC.append(0.0)

        if config.concParams.concProfile == "vonMises":
            ConcentrationA = self.VonMises(config.concParams.VonMisesVar, config.concParams.VonMisesMagnitude, config.concParams.VonMisesAOffset, self.length, config.simParams.locationStep)
            ConcentrationB = self.VonMises(config.concParams.VonMisesVar, config.concParams.VonMisesMagnitude, config.concParams.VonMisesBOffset, self.length, config.simParams.locationStep)
            ConcentrationC = self.VonMises(config.concParams.VonMisesVar, config.concParams.VonMisesMagnitude, config.concParams.VonMisesCOffset, self.length, config.simParams.locationStep)
        elif config.concParams.concProfile == "constant" or config.concParams.concProfile == "manaFromHeaven":
            ConcentrationA = self.constantConc(config.concParams.constantConc, self.length)
            ConcentrationB = self.constantConc(config.concParams.constantConc, self.length)
            ConcentrationC = self.constantConc(config.concParams.constantConc, self.length)
        elif config.concParams.concProfile == "constUpDown":
            ConcentrationA = self.constantUpDown(self.length, config.concParams.uphillLow, config.concParams.uphillHigh, config.concParams.locationA)
            ConcentrationB = self.constantUpDown(self.length, config.concParams.uphillLow, config.concParams.uphillHigh, config.concParams.locationB)
            ConcentrationC = self.constantUpDown(self.length, config.concParams.uphillLow, config.concParams.uphillHigh, config.concParams.locationC)

        CellLocations = config.cellMetaStats.cellLocations
        if len(CellLocations) == 0:
            CellLocations = range(self.length-1)
        Allcells = []
        for i in range(len(CellLocations)):
            currCell = cell.cell(newRand=1, Arec = config.cellStats.Arec, Brec = config.cellStats.Brec, Amol = config.cellStats.Amol, Bmol = config.cellStats.Bmol, config = config)
            Allcells.append(currCell)

        self.SimEnviornment = enviornment.enviornment(self.length, ConcentrationA, ConcentrationB, ConcentrationC, Allcells, CellLocations, newRand = 1, config = config)
        self.SimLength = config.simParams.simLength
        self.SimEnviornment.timeStep = config.simParams.simTimeStep
        return

    #If there was diffusion
    def runBeginningConc(self, time):
        for i in range(time):
            self.SimEnviornment.runConcentrationAdjusted()

    def moveRunReceptors(self):
        for i in range(self.SimLength):
            self.SimEnviornment.runCells()
            if i % 10 == 0:
                print(i)
                utils.plotAllCells(self.SimEnviornment, i)
        return self.SimEnviornment.Areceptors, self.SimEnviornment.Breceptors, self.SimEnviornment.allLocations

    def moveRunAB(self):
        Amol = []
        Bmol = []
        Acons =  []
        Bcons = []
        for i in range(self.SimLength):
            self.SimEnviornment.runCells()
            A,B, outA, outB = self.SimEnviornment.giveCellAB()
            self.SimEnviornment.decreaseSpace()

            Amol = Amol + A
            Bmol = Bmol + B
            Acons = Acons + outA
            Bcons = Bcons + outB

            #if i % 1 == 0:
                #print(i)
                #utils.plotAllCells(self.SimEnviornment, i)
        print(self.SimEnviornment.divisions)
        return Amol, Bmol, Acons, Bcons

    #records the movement of a cell in a run
    def moveRun(self):
        locs = []
        vel = 0.0
        multiple = 1.0
        totalCells = []
        divides_per_time = []
        realCells = []
        for i in range(self.SimLength):
            self.SimEnviornment.divides_per_time = 0
            self.SimEnviornment.runCells()
            print(self.SimEnviornment.time )
            locations = self.SimEnviornment.giveALlCellLocations()
            vels = self.SimEnviornment.giveAllCellVelocities()
            cellsAlive = self.SimEnviornment.cellsAlive()
            divides_per_time.append(self.SimEnviornment.divides_per_time*multiple)
            totalCells.append(cellsAlive * multiple)
            realCells.append(cellsAlive)

            if (self.SimEnviornment.decreaseSpace()):
                multiple = multiple * (cellsAlive / self.SimEnviornment.lowSpace)
            for j in range(len(vels)):
                vel = vel + vels[j]

            for j in range(len(locations)):
                locs.append(locations[j])
        print(vel)
        return locs, self.SimEnviornment.divides, self.SimEnviornment.splitLoc, self.SimEnviornment.Amol, self.SimEnviornment.Bmol, self.SimEnviornment.Aconcs, self.SimEnviornment.Bconcs, totalCells, divides_per_time, self.SimEnviornment.mol_times, realCells, self.SimEnviornment.Aconcentrations, self.SimEnviornment.Bconcentrations

    #simulates the enviornment in static concnetration and returns the division rate and MI of the environment
    def staicConcRunTime(self, printRun):
        retDivisions = []
        MIs = []
        Hx = []
        Hxgiveny = []
        ret2 = []
        entr = []
        entr_in = []
        roll = 10
        dataX1 = []
        dataX2 = []
        dataX3 = []
        dataX4 = []
        dataX5 = []
        dataX6 = []

        dataY1 = []
        dataY2 = []
        dataY3 = []
        dataY4 = []
        dataY5 = []
        dataY6 = []


        for i in range(roll):
            dataX1.append([])
            dataX2.append([])
            dataX3.append([])
            dataX4.append([])
            dataX5.append([])
            dataX6.append([])

            dataY1.append([])
            dataY2.append([])
            dataY3.append([])
            dataY4.append([])
            dataY5.append([])
            dataY6.append([])

        totalCells = []
        multiple = 1
        for i in range(self.SimLength/self.SimEnviornment.timeStep):
            #print(i)
            #print(self.SimEnviornment.cellsAlive())
            if i == 1:
                print(self.SimEnviornment.cellsAlive())
                #print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
            rollingEntropy = 0.0
            arr = self.SimEnviornment.Bconcentrations
            self.SimEnviornment.runCells()
            #print("made it here")
            if i % 100 == 0:
                #print(str(i) + " time")
                if printRun:
                    utils.plotAllCells(self.SimEnviornment, i)
                    #utils.plotReceptors(self.SimEnviornment, i)
            cellsAlive = self.SimEnviornment.cellsAlive()
            totalCells.append(cellsAlive*multiple)

            if(self.SimEnviornment.decreaseSpace()):
                multiple = multiple *(cellsAlive/self.SimEnviornment.lowSpace)
            #if i % 5 == 0:
            #    print(str(i) + " time")
            #    print(self.SimEnviornment.cellsAlive())
            if i % (self.SimLength-1) == 0 and i != 0:
                concs = self.SimEnviornment.giveAllCurrentConc()
                sum = 0.0
                for i in range(len(totalCells)):
                    sum += totalCells[i]/(len(totalCells))

                retDivisions.append(sum)

                dataX = []
                #print(self.SimEnviornment.cellsAlive())
                dataX1[i % roll] = concs[0]
                dataX2[i % roll] = concs[1]
                dataX3[i % roll] = concs[2]
                dataX4[i % roll] = concs[3]
                dataX5[i % roll] = concs[8]
                dataX6[i % roll] = concs[9]
                datx1 = []
                datx2 = []
                datx3 = []
                datx4 = []
                datx5 = []
                datx6 = []

                for j in range(len(dataX1)):
                    for k in range(len(dataX1[j])):
                        datx1.append(dataX1[j][k])
                    for k in range(len(dataX2[j])):
                        datx2.append(dataX2[j][k])
                    for k in range(len(dataX3[j])):
                        datx3.append(dataX3[j][k])
                    for k in range(len(dataX4[j])):
                        datx4.append(dataX4[j][k])
                    for k in range(len(dataX5[j])):
                        datx5.append(dataX5[j][k])
                    for k in range(len(dataX6[j])):
                        datx6.append(dataX6[j][k])

                dataX.append(datx1)
                dataX.append(datx2)
                dataX.append(datx3)
                dataX.append(datx4)
                #dataX.append(datx5)
                #dataX.append(datx6)

                dataY = []
                dataY1[i % roll] = concs[4]
                dataY2[i % roll] = concs[5]
                dataY3[i % roll] = concs[6]
                dataY4[i % roll] = concs[7]
                dataY5[i % roll] = concs[10]
                dataY6[i % roll] = concs[11]
                daty1 = []
                daty2 = []
                daty3 = []
                daty4 = []
                daty5 = []
                daty6 = []
                for j in range(len(dataY1)):
                    for k in range(len(dataY1[j])):
                        daty1.append(dataY1[j][k])
                    for k in range(len(dataY2[j])):
                        daty2.append(dataY2[j][k])
                    for k in range(len(dataY3[j])):
                        daty3.append(dataY3[j][k])
                    for k in range(len(dataY4[j])):
                        daty4.append(dataY4[j][k])
                    for k in range(len(dataY5[j])):
                        daty5.append(dataY5[j][k])
                    for k in range(len(dataY6[j])):
                        daty6.append(dataY6[j][k])
                dataY.append(daty1)
                dataY.append(daty2)
                dataY.append(daty3)
                dataY.append(daty4)
                #dataY.append(daty5)
                #dataY.append(daty6)

                testMI = MICalc.MICalc()
                #print(dataX)
                #print(dataY)
                ret = testMI.AltMI(dataX,dataY)
                ent_ret = testMI.entropy_allX(dataY)
                entr_in_ret = testMI.entropy_allX(dataX)
                #ret = 0;
                #
                print(ret)
                ret2.append(rollingEntropy)
                MIs.append(ret[0])
                Hx.append(ret[1])
                Hxgiveny.append(ret[2])
                entr.append(ent_ret)
                entr_in.append(entr_in_ret)
            if i % 1 == 0:
                self.SimEnviornment.resetDivisions()
        #print(MIs)
        #exit(-1)
        print(totalCells)
        move = self.SimEnviornment.giveMovement()
        sum1 = 0;
        for i in range(len(move)):
            sum1 = sum1 + math.fabs(move[i])
        sum1 = sum1/len(move)
        print(str(sum1) + "------------------------------------------------asdfasdgasfdgasdfvsadvsdavsdagvsfda------------------------------------------")
        print(MIs)
        print("MIs- -----------------------------------------------------------------------------------------------------------")
        return [totalCells, MIs, entr, entr_in, Hx, Hxgiveny]

    def staicConcRunTimeVar(self, printRun):
        retDivisions = []
        MIs = []
        Hx = []
        Hxgiveny = []
        ret2 = []
        entr = []
        entr_in = []
        roll = 10
        dataX1 = []
        dataX2 = []
        dataX3 = []
        dataX4 = []
        dataX5 = []
        dataX6 = []

        dataY1 = []
        dataY2 = []
        dataY3 = []
        dataY4 = []
        dataY5 = []
        dataY6 = []

        enviornment_Stats = []


        for i in range(roll):
            dataX1.append([])
            dataX2.append([])
            dataX3.append([])
            dataX4.append([])
            dataX5.append([])
            dataX6.append([])

            dataY1.append([])
            dataY2.append([])
            dataY3.append([])
            dataY4.append([])
            dataY5.append([])
            dataY6.append([])

        totalCells = []
        multiple = 1
        MI_vars = []
        for i in range(math.floor(self.SimLength/self.SimEnviornment.timeStep)):
            #print(i)
            #print(self.SimEnviornment.cellsAlive())
            #if i == 1:
                #print(self.SimEnviornment.cellsAlive())
                #print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
            rollingEntropy = 0.0
            arr = self.SimEnviornment.Bconcentrations
            self.SimEnviornment.runCells()
            #print("made it here")
            if i % 10 == 0:
                print(str(i) + " time")
                if printRun:
                    utils.plotAllCells(self.SimEnviornment, i)
                    utils.plotAllConcentrat(self.SimEnviornment, i)
            cellsAlive = self.SimEnviornment.cellsAlive()
            totalCells.append(cellsAlive*multiple)

            if(self.SimEnviornment.decreaseSpace()):
                multiple = multiple *(cellsAlive/self.SimEnviornment.lowSpace)
            #if i % 5 == 0:
            #    print(str(i) + " time")
            #    print(self.SimEnviornment.cellsAlive())
            if i == math.floor(self.SimLength/self.SimEnviornment.timeStep)-1:
                concs = self.SimEnviornment.giveAllCurrentConc()
                sum = 0.0
                for i in range(len(totalCells)):
                    sum += totalCells[i]/(len(totalCells))

                retDivisions.append(sum)

                dataX = []
                #print(self.SimEnviornment.cellsAlive())
                dataX1[i % roll] = concs[0]
                dataX2[i % roll] = concs[1]
                dataX3[i % roll] = concs[2]
                dataX4[i % roll] = concs[3]
                dataX5[i % roll] = concs[8]
                dataX6[i % roll] = concs[9]
                datx1 = []
                datx2 = []
                datx3 = []
                datx4 = []
                datx5 = []
                datx6 = []

                for j in range(len(dataX1)):
                    for k in range(len(dataX1[j])):
                        datx1.append(dataX1[j][k])
                    for k in range(len(dataX2[j])):
                        datx2.append(dataX2[j][k])
                    for k in range(len(dataX3[j])):
                        datx3.append(dataX3[j][k])
                    for k in range(len(dataX4[j])):
                        datx4.append(dataX4[j][k])
                    for k in range(len(dataX5[j])):
                        datx5.append(dataX5[j][k])
                    for k in range(len(dataX6[j])):
                        datx6.append(dataX6[j][k])

                dataX.append(datx1)
                dataX.append(datx2)
                dataX.append(datx3)
                dataX.append(datx4)
                #dataX.append(datx5)
                #dataX.append(datx6)

                dataY = []
                dataY1[i % roll] = concs[4]
                dataY2[i % roll] = concs[5]
                dataY3[i % roll] = concs[6]
                dataY4[i % roll] = concs[7]
                dataY5[i % roll] = concs[10]
                dataY6[i % roll] = concs[11]
                daty1 = []
                daty2 = []
                daty3 = []
                daty4 = []
                daty5 = []
                daty6 = []
                for j in range(len(dataY1)):
                    for k in range(len(dataY1[j])):
                        daty1.append(dataY1[j][k])
                    for k in range(len(dataY2[j])):
                        daty2.append(dataY2[j][k])
                    for k in range(len(dataY3[j])):
                        daty3.append(dataY3[j][k])
                    for k in range(len(dataY4[j])):
                        daty4.append(dataY4[j][k])
                    for k in range(len(dataY5[j])):
                        daty5.append(dataY5[j][k])
                    for k in range(len(dataY6[j])):
                        daty6.append(dataY6[j][k])
                dataY.append(daty1)
                dataY.append(daty2)
                dataY.append(daty3)
                dataY.append(daty4)
                #dataY.append(daty5)
                #dataY.append(daty6)

                testMI = MICalc.MICalc()
                #print(dataX)
                #print(dataY)
                print(len(self.SimEnviornment.giveAllCells()))
                ret = testMI.AltMI(dataX,dataY)
                ent_ret = testMI.entropy_allX(dataY)
                entr_in_ret = testMI.entropy_allX(dataX)
                #ret = 0;
                #
                ret2.append(rollingEntropy)
                MIs.append(ret[0])
                Hx.append(ret[1])
                Hxgiveny.append(ret[2])
                entr.append(ent_ret)
                entr_in.append(entr_in_ret)

            input_var = self.SimEnviornment.Recep_Var_Arr_input
            output_var = self.SimEnviornment.Recep_Var_Arr_Output

            if i == math.floor(self.SimLength/self.SimEnviornment.timeStep)-1:
                MIs_var_array = []
                MI_var = MICalc.MICalc()
                arr = []
                for k in range(len(input_var)):
                    if len(input_var[k][0]) == 0:
                        continue
                    input_var_MI = [input_var[k][0], input_var[k][1], input_var[k][2], input_var[k][3]]
                    output_var_MI = [output_var[k][0], output_var[k][1], output_var[k][2], output_var[k][3]]
                    var_length = len(input_var[k][0])
                    curr_MI = MI_var.AltMI(input_var_MI, output_var_MI)[0]
                    for j in range(var_length):
                        MIs_var_array.append(curr_MI)
                    arr.append(len(input_var[i][1]))
                mean = 0.0
                for k in range(len(MIs_var_array)):
                    mean = mean + MIs_var_array[k]
                mean = mean/len(MIs_var_array)
                final_var = 0.0
                for k in range(len(MIs_var_array)):
                    final_var = final_var + math.pow(math.fabs(MIs_var_array[k] - mean),2)
                final_var = math.sqrt(final_var/len(MIs_var_array))
                MI_vars.append(final_var)
                plt.grid()
                plt.hist(MIs_var_array, density = True, bins = 50)
                plt.clf()

            if i % 1 == 0:
                self.SimEnviornment.resetDivisions()

            enviornment_Stats.append(self.SimEnviornment.getAllStats())

        move = self.SimEnviornment.giveMovement()
        sum1 = 0;
        for i in range(len(move)):
            sum1 = sum1 + math.fabs(move[i])
        sum1 = sum1/len(move)
        return [totalCells, MIs, entr, entr_in, Hx, Hxgiveny, MI_vars, enviornment_Stats]

    #runs the simulation and calculates the MI, Division rate and calculated SNR
    def staicConcRunTimeSNR(self, printRun):
        retDivisions = []
        MIs = []
        SNRs = []
        for i in range(self.SimLength):
            #print(str(i) + " time")
            self.SimEnviornment.runCells()

            if i % 10 == 0:
                print(str(i) + " time")
                if printRun:
                    utils.plotAllCells(self.SimEnviornment, i)
                    utils.plotReceptors(self.SimEnviornment, i)
            if i % 1 == 0:
                concs = self.SimEnviornment.giveAllCurrentConc()
                allReceptors, concnentrations = self.SimEnviornment.giveAllReceptors()
                if self.SimEnviornment.cellsAlive() > 0:
                    retDivisions.append((self.SimEnviornment.giveDivisions()-self.SimEnviornment.giveDeaths())/self.SimEnviornment.cellsAlive())
                else:
                    retDivisions.append(0)
                dataX = []
                dataX.append(concs[0])
                dataX.append(concs[1])
                dataX.append(concs[2])
                dataX.append(concs[3])
                print(dataX)

                dataY = []
                dataY.append(concs[4])
                dataY.append(concs[5])
                dataY.append(concs[6])
                dataY.append(concs[7])
                testMI = MICalc.MICalc()
                ret = testMI.AltMI(dataX,dataY)
                MIs.append(ret)
                #SNR
                SNR = self.findSNRfromTable(allReceptors, concnentrations)
                SNRs.append(SNR)
        return [retDivisions, MIs, SNRs]

    def traditionalRun(self,time):
        self.runBeginningConc(time)
        frequency = self.concentrationFreq
        for i in range(self.SimLength):
            print("Time: " + str(i) + "-------------------------------------------------------------------------------")
            print(self.SimEnviornment.toString())
            self.SimEnviornment.runConcentrationAdjusted()
            self.SimEnviornment.runCells()
            if i % 10 == 0:
                utils.plotAllCells(self.SimEnviornment, i)
                utils.plotReceptors(self.SimEnviornment, i)
            if frequency > 0:
                if i % frequency == 0:
                    self.SimEnviornment.addAB(self.concentrationMag, self.locationA, self.locationB)

    def test(self, config, run, date):
        lam = config.simParams.lam
        gamma = config.simParams.gamma
        Aknoght = config.simParams.aknoght

        totalCells = []
        multiple = 1
        muadj = self.SimEnviornment.timeStep * lam

        cellLocations = []
        internalAB = []
        totalReceptors = []
        boundReceptors = []
        totalConcs = []
        cellMovement = []
        enviornmentConcs = []

        all_Aratios = []
        for i in range(10):
            all_Aratios.append([])

        for i in range(math.ceil(self.SimLength*(1/self.SimEnviornment.timeStep))):
            if i%((1/self.SimEnviornment.timeStep)*10) == 0:
                print("time: " + str(i*self.SimEnviornment.timeStep))
            self.SimEnviornment.runCells()

            if config.concParams.concProfile == "manaFromHeaven":
                # run mana from heaven
                if not config.simParams.presetBool:
                    rand1 = random.uniform(0, 1)
                    rand2 = random.uniform(0, 1)
                    if rand1 < muadj:
                        locationA = random.randint(0, self.SimEnviornment.length-1)
                        if config.simParams.presetSave:
                            config.simParams.presetA_time.append(i)
                            config.simParams.presetA_loc.append(locationA)
                        self.SimEnviornment.addA(Aknoght / self.SimEnviornment.locationStep, locationA)
                    if rand2 < muadj:
                        locationB = random.randint(0, self.SimEnviornment.length-1)
                        if config.simParams.presetSave:
                            config.simParams.presetB_time.append(i)
                            config.simParams.presetB_loc.append(locationB)
                        self.SimEnviornment.addB(Aknoght/self.SimEnviornment.locationStep, locationB)
                else:
                    if i in config.simParams.presetA_time:
                        loc = int(config.simParams.presetA_loc[config.simParams.presetA_time.index(i)])
                        self.SimEnviornment.addA(Aknoght / self.SimEnviornment.locationStep, loc)
                    if i in config.simParams.presetB_time:
                        loc = int(config.simParams.presetB_loc[config.simParams.presetB_time.index(i)])
                        self.SimEnviornment.addB(Aknoght / self.SimEnviornment.locationStep, loc)
                #run degradation
                degradeCoe = gamma * self.SimEnviornment.timeStep
                for j in range(len(self.SimEnviornment.Aconcentrations)):
                    self.SimEnviornment.Aconcentrations[j] = self.SimEnviornment.Aconcentrations[j]*(1-degradeCoe)
                    self.SimEnviornment.Bconcentrations[j] = self.SimEnviornment.Bconcentrations[j]*(1-degradeCoe)
                self.SimEnviornment.runConcentrationAdjusted()

            if config.outputFlags.concProfile:
                enviornmentConcs.append([])
                enviornmentConcs[i].append(self.SimEnviornment.Aconcentrations)
                enviornmentConcs[i].append(self.SimEnviornment.Bconcentrations)

            if config.outputFlags.cellLocations or config.outputFlags.internalAB or config.outputFlags.totalReceptors or config.outputFlags.boundReceptors or config.outputFlags.totalConcs or config.outputFlags.cellMovement:
                cells = self.SimEnviornment.positionHash
                if config.outputFlags.cellLocations:
                    cellLocations.append([])
                if config.outputFlags.internalAB:
                    internalAB.append([])
                if config.outputFlags.totalReceptors:
                    totalReceptors.append([])
                if config.outputFlags.boundReceptors:
                    boundReceptors.append([])
                if config.outputFlags.totalConcs:
                    totalConcs.append([])
                if config.outputFlags.cellMovement:
                    cellMovement.append([])
                for pos in cells:
                    cellsAtPos = cells[pos]
                    for cell in range(len(cellsAtPos)):
                        if config.outputFlags.cellLocations:
                            cellLocations[i].append(pos)
                        if config.outputFlags.internalAB:
                            internalAB[i].append([cellsAtPos[cell].Amol, cellsAtPos[cell].Bmol])
                        if config.outputFlags.totalReceptors:
                            totalReceptors[i].append([cellsAtPos[cell].Arec, cellsAtPos[cell].Brec])
                        if config.outputFlags.boundReceptors:
                            boundReceptors[i].append([cellsAtPos[cell].leftBoundArec, cellsAtPos[cell].rightBoundArec, cellsAtPos[cell].leftBoundBrec, cellsAtPos[cell].rightBoundBrec])
                        if config.outputFlags.totalConcs:
                            left = (pos-round(1/self.config.simParams.locationStep))% math.floor(self.length)
                            right = (pos+round(1/self.config.simParams.locationStep))% math.floor(self.length)
                            totalConcs[i].append([self.SimEnviornment.Aconcentrations[left], self.SimEnviornment.Aconcentrations[right] ,self.SimEnviornment.Bconcentrations[left], self.SimEnviornment.Bconcentrations[right]])
                        if config.outputFlags.cellMovement:
                            cellMovement[i].append(cellsAtPos[cell].vel)

            if (self.SimEnviornment.decreaseSpace()):
                if config.outputFlags.totalcells:
                    cellsAlive = self.SimEnviornment.cellsAlive()
                    totalCells.append(cellsAlive * multiple)
                    multiple = multiple * (cellsAlive / self.SimEnviornment.lowSpace)
            else:
                if config.outputFlags.totalcells:
                    cellsAlive = self.SimEnviornment.cellsAlive()
                    totalCells.append(cellsAlive * multiple)
            """""
            if not config.simParams.presetBool:
                if rand1 < muadj:
                    presetA.append(i)
                    locationA = random.randint(0, math.floor(self.SimEnviornment.length*(1/self.SimEnviornment.locationStep))-1)
                    presetA_loc.append(locationA)
                if rand2 < muadj:
                    presetB.append(i)
                    locationB = random.randint(0, math.floor(self.SimEnviornment.length*(1/self.SimEnviornment.locationStep))-1)
                    presetB_loc.append(locationB)
            else:
                for j in range(len(presetA)):
                    if i == presetA[j]:
                        self.SimEnviornment.addA(Aknoght / self.SimEnviornment.locationStep, presetA_loc[j])
                for j in range(len(presetB)):
                    if i == presetB[j]:
                        self.SimEnviornment.addB(Aknoght / self.SimEnviornment.locationStep, presetB_loc[j])
            """""

        totalCellsFile = None
        totalReceptorsFile = None
        boundReceptorsFile = None
        internalABFile = None
        totalConcsFile = None
        cellLocationsFile = None
        cellMovementFile = None
        enviornmentConcsFile = None
        if config.outputFlags.totalcells:
            totalCellsFile = utils.saveDataDate(run, date, "totalCells", totalCells, config.runOutputFlags.compressSave)
        if config.outputFlags.totalReceptors:
            totalReceptorsFile = utils.saveDataDate(run, date, "totalReceptors", totalReceptors, config.runOutputFlags.compressSave)
        if config.outputFlags.boundReceptors:
            boundReceptorsFile = utils.saveDataDate(run, date, "boundReceptors", boundReceptors, config.runOutputFlags.compressSave)
        if config.outputFlags.internalAB:
            internalABFile = utils.saveDataDate(run, date, "internalAB", internalAB, config.runOutputFlags.compressSave)
        if config.outputFlags.totalConcs:
            totalConcsFile = utils.saveDataDate(run, date, "totalConcs", totalConcs, config.runOutputFlags.compressSave)
        if config.outputFlags.cellLocations:
            cellLocationsFile = utils.saveDataDate(run, date, "cellLocations", cellLocations, config.runOutputFlags.compressSave)
        if config.outputFlags.cellMovement:
            cellMovementFile = utils.saveDataDate(run, date, "cellMovement", cellMovement, config.runOutputFlags.compressSave)
        if config.outputFlags.concProfile:
            enviornmentConcsFile = utils.saveDataDate(run, date, "enviornmentConcs", enviornmentConcs, config.runOutputFlags.compressSave)

        #print("Cells Alive At End: " + str(self.SimEnviornment.cellsAlive()))
        config.writeConfig("Data/" + date + "/" + run + "/config.cfg")
        out = output.output(runName = run,enviornmentConcsFile = enviornmentConcsFile, totalcellsFile = totalCellsFile, totalReceptorsFile = totalReceptorsFile, boundReceptorsFile = boundReceptorsFile, internalABFile = internalABFile, totalConcsFile = totalConcsFile, cellLocationsFile = cellLocationsFile, cellMovementFile = cellMovementFile)
        return out

    def indtest(self, config, run):
        lam = 0.4
        gamma = 0.2
        Aknoght = 400
        totalCells = []
        multiple = 1
        muadj = self.SimEnviornment.timeStep * lam

        cellLocations = []
        internalAB = []
        totalReceptors = []
        boundReceptors = []
        totalConcs = []
        cellMovement = []

        divisionsPerTime = []
        deathsPerTime = []

        all_Aratios = []
        for i in range(10):
            all_Aratios.append([])
        i = -1
        allDivisons = []
        allDeaths = []
        while self.SimEnviornment.cellsAlive() > 0:
            print("cells alive: " + str(self.SimEnviornment.cellsAlive()))
            i += 1
            if i%((1/self.SimEnviornment.timeStep)*1) == 0:
                print("time: " + str(i*self.SimEnviornment.timeStep))
                x = range(self.SimEnviornment.length)
                plt.bar(x, self.SimEnviornment.Aconcentrations, width=1)
                plt.bar(x, self.SimEnviornment.Bconcentrations, width=1)
                plt.hist(self.SimEnviornment.giveCellInArr(), 100)
                plt.show()

            self.SimEnviornment.runCells()
            allDivisons.append(self.SimEnviornment.divisions)
            allDeaths.append(self.SimEnviornment.giveDeaths())
            self.SimEnviornment.resetDivisions()

            if config.outputFlags.cellLocations or config.outputFlags.internalAB or config.outputFlags.totalReceptors or config.outputFlags.boundReceptors or config.outputFlags.totalConcs or config.outputFlags.cellMovement:
                cells = self.SimEnviornment.positionHash
                if config.outputFlags.cellLocations:
                    cellLocations.append([])
                if config.outputFlags.internalAB:
                    internalAB.append([])
                if config.outputFlags.totalReceptors:
                    totalReceptors.append([])
                if config.outputFlags.boundReceptors:
                    boundReceptors.append([])
                if config.outputFlags.totalConcs:
                    totalConcs.append([])
                if config.outputFlags.cellMovement:
                    cellMovement.append([])
                for pos in cells:
                    cellsAtPos = cells[pos]
                    for cell in range(len(cellsAtPos)):
                        if config.outputFlags.cellLocations:
                            cellLocations[i].append(pos)
                        if config.outputFlags.internalAB:
                            internalAB[i].append([cellsAtPos[cell].Amol, cellsAtPos[cell].Bmol])
                        if config.outputFlags.totalReceptors:
                            totalReceptors[i].append([cellsAtPos[cell].Arec, cellsAtPos[cell].Brec])
                        if config.outputFlags.boundReceptors:
                            boundReceptors[i].append([cellsAtPos[cell].leftBoundArec, cellsAtPos[cell].rightBoundArec, cellsAtPos[cell].leftBoundBrec, cellsAtPos[cell].rightBoundBrec])
                        if config.outputFlags.totalConcs:
                            left = (pos-round(1/self.config.simParams.locationStep))% math.floor(self.length)
                            right = (pos+round(1/self.config.simParams.locationStep))% math.floor(self.length)
                            totalConcs[i].append([self.SimEnviornment.Aconcentrations[left], self.SimEnviornment.Aconcentrations[right] ,self.SimEnviornment.Bconcentrations[left], self.SimEnviornment.Bconcentrations[right]])
                        if config.outputFlags.cellMovement:
                            cellMovement[i].append(cellsAtPos[cell].vel)

            #if (self.SimEnviornment.decreaseSpace()):
            #    if config.outputFlags.totalcells:
            #        cellsAlive = self.SimEnviornment.cellsAlive()
            #        totalCells.append(cellsAlive * multiple)
            #        multiple = multiple * (cellsAlive / self.SimEnviornment.lowSpace)
            #else:
            #    if config.outputFlags.totalcells:
            #        cellsAlive = self.SimEnviornment.cellsAlive()
            #        totalCells.append(cellsAlive * multiple)
            """""
            if not config.simParams.presetBool:
                if rand1 < muadj:
                    presetA.append(i)
                    locationA = random.randint(0, math.floor(self.SimEnviornment.length*(1/self.SimEnviornment.locationStep))-1)
                    presetA_loc.append(locationA)
                if rand2 < muadj:
                    presetB.append(i)
                    locationB = random.randint(0, math.floor(self.SimEnviornment.length*(1/self.SimEnviornment.locationStep))-1)
                    presetB_loc.append(locationB)
            else:
                for j in range(len(presetA)):
                    if i == presetA[j]:
                        self.SimEnviornment.addA(Aknoght / self.SimEnviornment.locationStep, presetA_loc[j])
                for j in range(len(presetB)):
                    if i == presetB[j]:
                        self.SimEnviornment.addB(Aknoght / self.SimEnviornment.locationStep, presetB_loc[j])
            """""

        totalCellsFile = None
        totalReceptorsFile = None
        boundReceptorsFile = None
        internalABFile = None
        totalConcsFile = None
        cellLocationsFile = None
        cellMovementFile = None
        if config.outputFlags.totalcells:
            totalCellsFile = utils.saveDataDate(run, "totalCells_" + run, totalCells, config.runOutputFlags.compressSave)
        if config.outputFlags.totalReceptors:
            totalReceptorsFile = utils.saveDataDate(run, "totalReceptors_" + run, totalReceptors, config.runOutputFlags.compressSave)
        if config.outputFlags.boundReceptors:
            boundReceptorsFile = utils.saveDataDate(run, "boundReceptors_" + run, boundReceptors, config.runOutputFlags.compressSave)
        if config.outputFlags.internalAB:
            internalABFile = utils.saveDataDate(run, "internalAB_" + run, internalAB, config.runOutputFlags.compressSave)
        if config.outputFlags.totalConcs:
            totalConcsFile = utils.saveDataDate(run, "totalConcs_" + run, totalConcs, config.runOutputFlags.compressSave)
        if config.outputFlags.cellLocations:
            cellLocationsFile = utils.saveDataDate(run, "cellLocations_" + run, cellLocations, config.runOutputFlags.compressSave)
        if config.outputFlags.cellMovement:
            cellMovementFile = utils.saveDataDate(run, "cellMovement_" + run, cellMovement, config.runOutputFlags.compressSave)
        print("Cells Alive At End: " + str(self.SimEnviornment.cellsAlive()))
        out = output.output(totalcellsFile = totalCellsFile, totalReceptorsFile = totalReceptorsFile, boundReceptorsFile = boundReceptorsFile, internalABFile = internalABFile, totalConcsFile = totalConcsFile, cellLocationsFile = cellLocationsFile, cellMovementFile = cellMovementFile)
        return out, allDivisons, allDeaths

