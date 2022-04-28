import cell
import utils
import math
import copy
import random

class enviornment(object):
    #enviornment parameters: length, cells, cell positions, concentrations
    #cell positions define the left of the cell, the right is the position + 1, cells can only go to length - 2
    Aabsorbed = 0
    Babsorbed = 0

    distancetraveled = 0

    space = 10000
    lowSpace = 9000

    def findCurrentCellCount(self):
        allCells = self.positionHash
        sum = 0
        for i in allCells:
            sum += len(allCells[i])
        return sum

    def __init__(self, length, Aconcentrations, Bconcentrations, DiffusionCoe, setOfCells, cellPositions, largestID, newRand, Cconcentrations):
        self.positionHash = utils.createPositionHash(length)
        for i in range(len(cellPositions)):
            self.positionHash[cellPositions[i]].append(setOfCells[i])
        self.time = 0
        self.Aconcentrations = Aconcentrations
        print(len(Aconcentrations))
        self.Bconcentrations = Bconcentrations
        self.Cconcentrations = Cconcentrations

        self.Adiffusion = DiffusionCoe #1x^2/time = 1^2/1 = 1
        self.Bdiffusion = DiffusionCoe #1x^2/time = 1^2/1 = 1
        self.NumberOfCells = len(setOfCells)
        self.length = length
        self.CellDeaths = 0
        self.currentID = largestID + 1
        self.totalCellCount = length
        self.totalDivisions = 0
        self.alphaA = utils.highest(Aconcentrations)
        self.alphaB = utils.highest(Bconcentrations)

        self.allGradientsAl = []
        self.allGradientsBl = []
        self.allGradientsAr = []
        self.allGradientsBr = []

        self.allBoundAl = []
        self.allBoundBl = []
        self.allBoundAr = []
        self.allBoundBr = []
        self.cellVelocities = []

        self.fullDivide = True

        self.Amol = []
        self.Bmol = []
        #self.time = -1
        self.mol_times = []
        self.Aconcs = []
        self.Bconcs = []

        self.divisions = 0
        self.deaths = 0

        self.logdoublingtimes = []

        self.allLocations = []

        self.pastCellCount = self.findCurrentCellCount()

        self.newRand = newRand

        self.fullDie = True

        self.divideLoc = []

        self.splitLoc = []

        self.Areceptors = []
        self.Breceptors = []

        self.IntA = []
        self.IntB = []

        self.divides = []
        self.divides_per_time = 0

        self.Recep_Var_Arr_input = []
        self.Recep_Var_Arr_Output = []

    def addlogdoubletime(self):
        currentCellCount = self.findCurrentCellCount()
        difference = currentCellCount - self.pastCellCount
        self.pastCellCount = currentCellCount
        if difference != 0.0:
            self.logdoublingtimes.append(currentCellCount/difference)

    def addAB(self, magnitude, locA, locB):
        if locA < 0:
            locA = self.newRand.getNewInt(0, self.length)
        if locB < 0:
            locB = self.newRand.getNewInt(0, self.length)
        self.Aconcentrations[locA] += magnitude
        self.Bconcentrations[locB] += magnitude

    def absorbAB(self, position, A, B):
        if A > 0.0:
            A1 = self.Aconcentrations[position]
            A2 = self.Aconcentrations[position + 1]
            newA1 = utils.checkzero(A1 - ((A1 / (A1 + A2)) * (A / 1000.0)))
            newA2 = utils.checkzero(A2 - ((A2 / (A1 + A2)) * (A / 1000.0)))
            self.Aconcentrations[position] = newA1
            self.Aconcentrations[position + 1] = newA2

        if B > 0.0:
            B1 = self.Bconcentrations[position]
            B2 = self.Bconcentrations[position + 1]
            newB1 = utils.checkzero(B1 - ((B1 / (B1 + B2)) * (B / 1000.0)))
            newB2 = utils.checkzero(B2 - ((B2 / (B1 + B2)) * (B / 1000.0)))
            self.Bconcentrations[position] = newB1
            self.Bconcentrations[position + 1] = newB2

    def runCell(self, position, cell):

        self.Amol.append(cell.Amol + cell.ATP)
        self.Bmol.append(cell.Bmol + cell.ATP)
        self.mol_times.append(self.time)
        self.divides.append(0)
        self.Aconcs.append(self.Aconcentrations[position])
        self.Bconcs.append(self.Bconcentrations[position])
        left = (position - 1) % self.length
        right = (position + 1) % self.length
        self.allGradientsAl.append(self.Aconcentrations[left])
        self.allGradientsBl.append(self.Bconcentrations[left])
        self.allGradientsAr.append(self.Aconcentrations[right])
        self.allGradientsBr.append(self.Bconcentrations[right])
        self.IntA.append(cell.Amol + cell.ATP)
        self.IntB.append(cell.Bmol + cell.ATP)

        self.divideLoc.append(0)

        [absorbA, absorbB] = cell.Absorb(self.Aconcentrations[left], self.Aconcentrations[right], self.Bconcentrations[left], self.Bconcentrations[right])
        self.Aabsorbed += absorbA
        self.Babsorbed += absorbB

        velocity = cell.velocity(self.Aconcentrations[left], self.Aconcentrations[right], self.Bconcentrations[left], self.Bconcentrations[right], self.alphaA, self.alphaB)

        self.allLocations.append(position)

        if position + velocity < 0:
            self.splitLoc.append(-1)
        elif position + velocity > self.length-1:
            self.splitLoc.append(1)
        else:
            self.splitLoc.append(0)

        self.cellVelocities.append(velocity)
        self.allBoundAl.append(cell.leftBoundArec)
        self.allBoundBl.append(cell.leftBoundBrec)
        self.allBoundAr.append(cell.rightBoundArec)
        self.allBoundBr.append(cell.rightBoundBrec)

        self.IntA.append(cell.Amol + cell.ATP)
        self.IntB.append(cell.Bmol + cell.ATP)

        #do this only after calculating the velocities
        self.Areceptors.append(cell.Arec)
        self.Breceptors.append(cell.Brec)
        cell.changeReceptors()

        if cell.consumption(velocity, self.fullDie) == False:
            self.divideLoc.append(0)
            self.Amol.append(cell.Amol + cell.ATP)
            self.Bmol.append(cell.Bmol + cell.ATP)
            self.mol_times.append(self.time)
            self.divides.append(0)
            self.Aconcs.append(self.Aconcentrations[position])
            self.Bconcs.append(self.Bconcentrations[position])
            return [0, -1]
        cell.incrementDis(math.fabs(velocity))
        if cell.divide() == True:
            self.divideLoc.append(1)
            self.Amol.append(cell.Amol + cell.ATP)
            self.Bmol.append(cell.Bmol + cell.ATP)
            self.mol_times.append(self.time)
            self.divides.append(1)
            self.Aconcs.append(self.Aconcentrations[position])
            self.Bconcs.append(self.Bconcentrations[position])
            return [velocity, 1]
        self.divideLoc.append(0)
        #.append(cell.Amol + cell.ATP)
        #self.Bmol.append(cell.Bmol + cell.ATP)
        self.Amol.append(cell.Amol + cell.ATP)
        self.Bmol.append(cell.Bmol + cell.ATP)
        self.divides.append(0)
        self.Aconcs.append(self.Aconcentrations[position])
        self.Bconcs.append(self.Bconcentrations[position])
        self.mol_times.append(self.time)
        return [velocity, 0]

    def changeLen(self, arr, newLen):
        newArr = []
        for i in range(0, newLen):
            newArr.append(arr[i])
        return newArr

    def decreaseSpace(self):
        AllCells = self.positionHash
        cellCount = 0.0
        for cells in AllCells:
            cellCount += float(len(AllCells[cells]))
        if cellCount < self.space:
            return False
        prob = self.lowSpace/cellCount
        for cells in AllCells:
            newLen = math.floor(len(AllCells[cells])*prob)
            newArr = self.changeLen(AllCells[cells], newLen)
            AllCells[cells] = newArr
        newLen = 0
        for cells in AllCells:
            newLen += len(AllCells[cells])
        return True

    def runCells(self):
        #self.time += 1
        newPositionHash = utils.createPositionHash(self.length)
        for pos in self.positionHash:
            cellslen = len(self.positionHash[pos])
            for i in range(cellslen):
                #print("another cell")
                movingCell = self.positionHash[pos][i]
                arr = self.runCell(pos, movingCell)
                #print("made it here")
                #it died
                if arr[1] == -1:
                    self.deaths += 1
                    continue
                newPosition = (pos + arr[0]) % self.length
                movingCell.incrementDis(math.fabs(newPosition - pos))
                #it didn't divide
                if arr[1] == 0:
                    #self.divides.append(0)
                    newPositionHash[newPosition].append(movingCell)
                #it divided
                if arr[1] == 1:
                    #self.divides.append(1)
                    self.divides_per_time += 1
                    self.divisions += 1
                    stats = movingCell.divideStats()
                    #self.Amol.append(stats[3] + stats[5])
                    #self.Bmol.append(stats[4] + stats[5])
                    Arec = stats[0]
                    Brec = stats[1]
                    MaxRec = stats[2]
                    Amol = stats[3]
                    Bmol = stats[4]
                    ATP = stats[5]
                    biomass = stats[6]
                    generation = stats[7]
                    distancetraveled = stats[8]

                    cell1 = cell.cell(Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distancetraveled, self.currentID, self.newRand)
                    cell1.AbsorbtionRate = movingCell.AbsorbtionRate
                    cell1.ReceptorConsumptionRate = movingCell.ReceptorConsumptionRate
                    cell1.survivalCost = movingCell.survivalCost
                    cell1.VelocityMultiplier = movingCell.VelocityMultiplier
                    cell1.mutate = movingCell.mutate
                    cell1.decisiontype = movingCell.decisiontype
                    cell1.noise = movingCell.noise
                    cell1.receptor_mode = movingCell.receptor_mode
                    newPositionHash[newPosition].append(cell1)

                    if self.fullDivide:
                        cell2 = cell.cell(Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distancetraveled, self.currentID + 1, self.newRand)
                        cell2.AbsorbtionRate = movingCell.AbsorbtionRate
                        cell2.ReceptorConsumptionRate = movingCell.ReceptorConsumptionRate
                        cell2.survivalCost = movingCell.survivalCost
                        cell2.VelocityMultiplier = movingCell.VelocityMultiplier
                        cell2.mutate = movingCell.mutate
                        cell2.decisiontype = movingCell.decisiontype
                        cell2.noise = movingCell.noise
                        cell2.receptor_mode = movingCell.receptor_mode
                        newPositionHash[newPosition].append(cell2)

                    self.totalCellCount += 1
                    self.totalDivisions += 1
                    self.currentID += 2
        Recep_Var_Arr_input = []
        Recep_Var_Arr_Output = []
        for i in range(401):
            Recep_Var_Arr_input.append([])
            Recep_Var_Arr_Output.append([])

            Recep_Var_Arr_input[i].append([])
            Recep_Var_Arr_input[i].append([])
            Recep_Var_Arr_input[i].append([])
            Recep_Var_Arr_input[i].append([])

            Recep_Var_Arr_Output[i].append([])
            Recep_Var_Arr_Output[i].append([])
            Recep_Var_Arr_Output[i].append([])
            Recep_Var_Arr_Output[i].append([])

        for pos in self.positionHash:
            cellslen = len(self.positionHash[pos])
            for i in range(cellslen):
                currCell = self.positionHash[pos][i]
                Arec = currCell.Arec
                #print(Arec)
                #print(len(Recep_Var_Arr_input))
                #print(len(Recep_Var_Arr_input[Arec]))
                Recep_Var_Arr_input[Arec][0].append(currCell.AconLeft)
                Recep_Var_Arr_input[Arec][1].append(currCell.AconRight)
                Recep_Var_Arr_input[Arec][2].append(currCell.BconLeft)
                Recep_Var_Arr_input[Arec][3].append(currCell.BconRight)

                Recep_Var_Arr_Output[Arec][0].append(currCell.leftBoundArec)
                Recep_Var_Arr_Output[Arec][1].append(currCell.rightBoundArec)
                Recep_Var_Arr_Output[Arec][2].append(currCell.leftBoundBrec)
                Recep_Var_Arr_Output[Arec][3].append(currCell.rightBoundBrec)

        self.Recep_Var_Arr_input = Recep_Var_Arr_input
        self.Recep_Var_Arr_Output = Recep_Var_Arr_Output

        self.positionHash = newPositionHash
        self.time += 1
        self.addlogdoubletime()
        #self.decreaseSpace()

    def runConcentration(self, difusionCo):
        newAcon = []
        newBcon = []
        for i in range(self.length):
            if i == 0:
                newAcon.append(self.Aconcentrations[i] + (difusionCo *
                                                        (self.Aconcentrations[i + 1] - self.Aconcentrations[i])))
                newBcon.append(self.Bconcentrations[i] + (difusionCo *
                                                        (self.Bconcentrations[i + 1] - self.Bconcentrations[i])))
                continue
            elif i == self.length - 1:
                newAcon.append(self.Aconcentrations[i] + (difusionCo *
                                                          (self.Aconcentrations[i - 1] - self.Aconcentrations[i])))
                newBcon.append(self.Bconcentrations[i] + (difusionCo *
                                                          (self.Bconcentrations[i - 1] - self.Bconcentrations[i])))
                continue
            newAcon.append(self.Aconcentrations[i] + (difusionCo *
                                                      (self.Aconcentrations[i+1] - 2*self.Aconcentrations[i] + self.Aconcentrations[i-1])))
            newBcon.append(self.Bconcentrations[i] + (difusionCo *
                                                      (self.Bconcentrations[i+1] - 2*self.Bconcentrations[i] + self.Bconcentrations[i-1])))
        self.Aconcentrations = newAcon
        self.Bconcentrations = newBcon

    def runConcentrationAdjusted(self):
        TimesRun = math.ceil(self.Adiffusion/0.5)
        adjustedAdiffusion = 1.0/TimesRun
        for i in range(TimesRun):
            self.runConcentration(adjustedAdiffusion)

    def toString(self):
        ret = ""
        for i in range(self.length):
            ret += ", " + str(len(self.positionHash[i]))
        return ret

    def AconcnentratToString(self):
        ret = ""
        for i in range(len(self.Aconcentrations)):
            ret += ", " + str(self.Aconcentrations[i])
        return ret

    def toArr(self):
        arr = []
        for i in range(self.length):
            arr.append(len(self.positionHash[i]))
        return arr

    def returnEfficency(self):
        AB = float(self.Aabsorbed + self.Babsorbed)
        d = float(self.divisions)
        return d/AB

    def returnAverageDistanceTraveled(self):
        j = 0.0
        averageDistance = 0.0
        for pos in self.positionHash:
            list = self.positionHash[pos]
            length = len(list)
            totaldistance = 0.0
            for i in list:
                totaldistance += float(i.distanceTrav)

            if len(list) > 0:
                averageDistance += totaldistance / length
                j+= 1.0
        if j == 0:
            j = 1
        return averageDistance/j

    def returnAverageRec(self):
        finalAverage = []
        for pos in self.positionHash:
            listA = self.positionHash[pos]
            average = [0.0, 0.0]
            if len(listA) == 0:
                finalAverage.append(average)
                continue
            AverageA = 0.0
            AverageB = 0.0
            for i in listA:
                AverageA += i.Arec
                AverageB += i.Brec
            AverageA = AverageA/len(listA)
            AverageB = AverageB/len(listA)
            average = [AverageA,AverageB]
            finalAverage.append(average)
        return finalAverage

    def cellsAlive(self):
        ret = 0
        for i in self.positionHash:
            ret += len(self.positionHash[i])
        return ret

    def returnAverageDoublingTime(self):
        sum = 0.0
        for i in range(len(self.logdoublingtimes)):
            sum += self.logdoublingtimes[i]
        return sum/len(self.logdoublingtimes)

    def giveAllConc(self):
        XAl = self.allGradientsAl
        XBl = self.allGradientsBl
        XAr = self.allGradientsAr
        XBr = self.allGradientsBr
        YAl = self.allBoundAl
        YBl = self.allBoundBl
        YAr = self.allBoundAr
        YBr = self.allBoundBr
        return XAl, XBl, XAr, XBr, YAl, YBl, YAr, YBr

    def giveAllCurrentConc(self):
        XAl = []
        XBl = []
        XAr = []
        XBr = []
        XCl = []
        XCr = []

        YAl = []
        YBl = []
        YAr = []
        YBr = []
        YCl = []
        YCr = []
        for i in range(self.length):
            Cells = self.positionHash[i]
            for j in Cells:
                YAl.append(j.leftBoundArec)
                YBl.append(j.leftBoundBrec)
                YAr.append(j.rightBoundArec)
                YBr.append(j.rightBoundBrec)
                #YCl.append(utils.calculateBoundKinetic(100,self.Cconcentrations[(i-1) %self.length], 2, 1))
                #YCr.append(utils.calculateBoundKinetic(100,self.Cconcentrations[(i+1) %self.length], 2, 1))

                XAl.append(self.Aconcentrations[(i-1) % self.length])
                XBl.append(self.Bconcentrations[(i + 1) % self.length])
                XAr.append(self.Aconcentrations[(i - 1) % self.length])
                XBr.append(self.Bconcentrations[(i + 1) % self.length])
                XCl.append(self.Cconcentrations[(i-1) % self.length])
                XCr.append(self.Cconcentrations[(i+1) % self.length])
        return XAl, XBl, XAr, XBr, YAl, YBl, YAr, YBr, XCl, XCr, YCl, YCr

    def giveAllReceptors(self):
        recs = []
        concs = []
        for i in range(self.length):
            cells = self.positionHash[i]
            for j in cells:
                recs.append(int(j.Arec/2))
                concs.append(self.Aconcentrations[(i-1) % self.length])
                recs.append(int(j.Arec / 2))
                concs.append(self.Aconcentrations[(i+1) % self.length])
                recs.append(int(j.Brec / 2))
                concs.append(self.Bconcentrations[(i-1) % self.length])
                recs.append(int(j.Brec / 2))
                concs.append(self.Bconcentrations[(i+1) % self.length])
        return recs, concs

    def giveDivisions(self):
        ret = self.divisions
        return ret

    def resetDivisions(self):
        self.divisions = 0

    def giveDeaths(self):
        ret = self.deaths
        self.deaths = 0
        return ret

    def resetConc(self):
        self.allGradientsAl = []
        self.allGradientsBl = []
        self.allGradientsAr = []
        self.allGradientsBr = []
        self.allBoundAl = []
        self.allBoundBl = []
        self.allBoundAr = []
        self.allBoundBr = []
        return

    def giveMovement(self):
        return self.cellVelocities

    def giveALlCellLocations(self):
        ret = []
        for i in range(self.length):
            cells = self.positionHash[i]
            for j in cells:
                ret.append(i)
        return ret

    def giveAllCellVelocities(self):
        ret = []
        for i in range(self.length):
            cells = self.positionHash[i]
            for j in cells:
                ret.append(j.vel)
        return ret

    def giveAllCells(self):
        ret = []
        for i in range(self.length):
            cells = self.positionHash[i]
            for j in cells:
                ret.append(j)
        return ret

    def giveCellAB(self):
        A = 0
        B = 0
        outA = 0.0
        outB = 0.0
        Amol1 = []
        Bmol1 = []
        Aconc = []
        Bconc = []
        for pos in self.positionHash:
            for i in range(len(self.positionHash[pos])):
                cell = self.positionHash[pos][i]
                cell.Amol + cell.ATP, cell.Bmol + cell.ATP, self.Aconcentrations[pos], self.Bconcentrations[pos]
                Amol1.append(cell.Amol+cell.ATP)
                Bmol1.append(cell.Bmol+cell.ATP)
                Aconc.append(self.Aconcentrations[pos])
                Bconc.append(self.Bconcentrations[pos])
        return Amol1, Bmol1, Aconc, Bconc