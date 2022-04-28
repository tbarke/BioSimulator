import cell
import enviornment
import MICalc
import utils
import random
import math
from scipy import special
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

    def VonMises(self, var, maximum, offset, length):
        length = length
        arr = []
        kappa = 1/var
        mu = 2 * math.pi / length * offset
        bessel = (special.iv(0, kappa))
        sum = 0.0
        increment = 2 * math.pi / length
        for i in range(length):
            index = i - length/2
            x = index * increment
            numerator = math.exp(kappa * math.cos(x - mu))
            arr.append(numerator / (2 * math.pi * bessel))
            sum += numerator / (2 * math.pi * bessel)
        arr = (arr / sum)*maximum
        return arr

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
    def __init__(self, ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams, SimParams, newRand):
        length = SimParams[0]
        #concentrations
        diffusionCoefficent = ConcParams[0]
        absorb = False

        if ConcParams[1] == "absorb":
            absorb = True
        self.concentrationFreq = ConcParams[2]
        self.concentrationMag = ConcParams[3]
        self.locationA = ConcParams[4]
        self.locationB = ConcParams[5]
        self.locationC = ConcParams[12]
        ConcentrationA = []
        ConcentrationB = []
        ConcentrationC = []
        for i in range(length):
            ConcentrationA.append(0.0)
            ConcentrationB.append(0.0)
            ConcentrationC.append(0.0)
        if self.locationA > -1:
            ConcentrationA[self.locationA] = self.concentrationMag
            ConcentrationB[self.locationB] = self.concentrationMag
            ConcentrationC[self.locationC] = self.concentrationMag
        else:
            ConcentrationA[newRand.getNewInt(0, length)] = self.concentrationMag
            ConcentrationB[newRand.getNewInt(0, length)] = self.concentrationMag
            ConcentrationC[newRand.getNewInt(0, length)] = self.concentrationMag
        if ConcParams[6] == "static":
            ConcentrationA = self.VonMises(ConcParams[7], ConcParams[8], ConcParams[9], length)
            ConcentrationB = self.VonMises(ConcParams[7], ConcParams[8], ConcParams[10], length)
            ret = "concA = ["
            ret2 = "concB = ["
            for i in range(0, len(ConcentrationA)-1):
                ret = ret + str(ConcentrationA[i]) + ", "
                ret2 = ret2 + str(ConcentrationB[i]) + ", "
            ret = ret + str(ConcentrationA[len(ConcentrationA)-1]) + "]"
            ret2 = ret2 + str(ConcentrationB[len(ConcentrationB)-1]) + "]"
            print(ret)
            print(ret2)
            #print(ret)
            #exit(-1)
            #exit(-1)
            ConcentrationC = self.VonMises(ConcParams[7], ConcParams[8], ConcParams[11], length)

        #cell locations data
        CellLocations = CellLocationStats[0]
        if len(CellLocations) == 0:
            CellLocations = range(length-1)
        largestID = CellStats[9]
        Allcells = []
        for i in range(len(CellLocations)):
            currCell = cell.cell(CellStats[0], CellStats[1], CellStats[2], CellStats[3], CellStats[4], CellStats[5], CellStats[6], CellStats[7], CellStats[8], largestID, newRand)
            currCell.AbsorbtionRate = CellMetaStats[0]
            currCell.ReceptorConsumptionRate = CellMetaStats[1]
            currCell.survivalCost = CellMetaStats[2]
            currCell.VelocityMultiplier = CellMetaStats[3]
            currCell.noise = CellMetaStats[6]
            currCell.receptor_mode = CellMetaStats[7]
            currCell.VelocityMultiplier = CellMetaStats[8]
            if CellMetaStats[4] == "non":
                currCell.mutate = False
            if CellMetaStats[5] == "non":
                currCell.decisiontype = "non"
            Allcells.append(currCell)
            largestID += 1

        self.SimEnviornment = enviornment.enviornment(length, ConcentrationA, ConcentrationB, diffusionCoefficent, Allcells, CellLocations, largestID, newRand, ConcentrationC)
        self.SimEnviornment.fullDivide = EnviornmentParams[1]
        self.SimEnviornment.fullDie = EnviornmentParams[2]
        self.SimLength = EnviornmentParams[0]
        self.dataRecord = SimParams
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
        for i in range(self.SimLength):
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
        for i in range(self.SimLength):
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

            input_var = self.SimEnviornment.Recep_Var_Arr_input
            output_var = self.SimEnviornment.Recep_Var_Arr_Output

            if i % 100 == 0:
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
                    #print(len(input_var[i][1]))
                mean = 0.0
                for k in range(len(MIs_var_array)):
                    mean = mean + MIs_var_array[k]
                mean = mean/len(MIs_var_array)
                final_var = 0.0
                for k in range(len(MIs_var_array)):
                    final_var = final_var + math.pow(math.fabs(MIs_var_array[k] - mean),2)
                final_var = math.sqrt(final_var/len(MIs_var_array))
                MI_vars.append(final_var)
                #print("here")


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
        #print(MI_vars)
        #exit(-1)
        return [totalCells, MIs, entr, entr_in, Hx, Hxgiveny, MI_vars]



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