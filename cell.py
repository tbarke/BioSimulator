import math
import utils
import numpy as np

class cell(object):
#cell attributes: A receptors, B receptors, A molecules, B molecules, Bound A receptors, and Bound B receptors
    biomassDivide = 10
    AbsorbtionRate = 1.0
    ReceptorConsumptionRate = 0.10
    survivalCost = 2 #ATP
    VelocityMultiplier = 1.0
    mutate = True
    decisiontype = "naive"
    noise = 0
    receptor_mode = "simulate"
    def __init__(self, Arec, Brec, MaxReceptors, Amol, Bmol, ATP, biomass, generation, distanceTrav, ID, newRand):
        self.Arec = Arec
        self.Brec = Brec
        self.MaxReceptors = MaxReceptors
        self.Amol = Amol
        self.Bmol = Bmol
        self.ATP = ATP
        self.generation = generation
        self.distanceTrav = distanceTrav
        self.biomass = biomass
        self.ID = ID

        #bound receptors, get rewritten every time step
        self.leftBoundArec = 0
        self.rightBoundArec = 0
        self.leftBoundBrec = 0
        self.rightBoundBrec = 0

        self.AconLeft = 0
        self.AconRight = 0
        self.BconLeft = 0
        self.BconRight = 0

        self.newRand = newRand
        self.vel = 0


    def Absorb(self, AconLeft, AconRight, BconLeft, BconRight):
        retA = utils.calculateAbsorbtion((AconLeft + AconRight) / 2, self.AbsorbtionRate)
        retB = utils.calculateAbsorbtion((BconLeft + BconRight) / 2, self.AbsorbtionRate)
        self.Amol = self.Amol + retA
        self.Bmol = self.Bmol + retB
        complex = utils.ATPcomplex(self.Amol, self.Bmol, 1, 1)
        self.ATP += complex[0]
        self.Amol -= complex[1]
        self.Bmol -= complex[2]
        return [retA, retB]

    #returns the velocity of the cell
    def velocity(self, AconLeft, AconRight, BconLeft, BconRight, alphaA, alphaB):
        leftArec = math.floor(self.Arec/2)
        rightArec = leftArec
        leftBrec = math.floor(self.Brec/2)
        rightBrec = leftBrec
        dissocociation = 2
        self.AconLeft = AconLeft
        self.AconRight = AconRight
        self.BconLeft = BconLeft
        self.BconRight = BconRight
        self.leftBoundArec = utils.calculateBoundKinetic(leftArec, AconLeft, dissocociation, self.newRand, self.noise, self.receptor_mode)
        self.rightBoundArec = utils.calculateBoundKinetic(rightArec, AconRight, dissocociation, self.newRand, self.noise, self.receptor_mode)
        self.leftBoundBrec = utils.calculateBoundKinetic(leftBrec, BconLeft, dissocociation, self.newRand, self.noise, self.receptor_mode)
        self.rightBoundBrec = utils.calculateBoundKinetic(rightBrec, BconRight, dissocociation, self.newRand, self.noise, self.receptor_mode)

        #print(self.leftBoundArec)
        #print(self.leftBoundBrec)
        #print(self.rightBoundArec)
        #print(self.rightBoundBrec)

        #adopting right positive notation
        vel = math.floor(((self.rightBoundArec + self.rightBoundBrec) - (
                    self.leftBoundArec + self.leftBoundBrec)) * self.VelocityMultiplier)
        #vel = math.floor(((self.leftBoundArec + self.leftBoundBrec) - (self.rightBoundArec + self.rightBoundBrec)) * self.VelocityMultiplier)
        #print(vel)
        #exit(-1)
        self.vel = vel
        rand = np.random.randint(-50,50, 1)[0]
        #print(rand)
        if vel > 10:
            vel = 10
        elif vel < -10:
            vel = -10
        return rand #math.floor(((self.rightBoundArec + self.rightBoundBrec) - (self.leftBoundArec + self.leftBoundBrec))*self.VelocityMultiplier)

    def consumption(self, velocity, fullDie):
        AllCost = self.survivalCost
        self.ATP -= AllCost
        if self.ATP < 0:
            self.Amol += self.ATP
            self.Bmol += self.ATP
            if self.Amol < 0:
                self.Amol = 0
            if self.Bmol < 0:
                self.Bmol = 0
            self.ATP = 0
            if fullDie:
                #print("this happened")
                return False
        return True

    def incrementDis(self, distance):
        self.distanceTrav += distance

    def divide(self):
        self.biomass += 1
        if self.ATP > self.survivalCost*5:
            return True
        return False

    def changeReceptors(self):
        #self.decisiontype = "naive"
        interA = self.Amol + self.ATP
        interB = self.Bmol + self.ATP
        exterA = self.leftBoundArec + self.rightBoundArec
        exterB = self.leftBoundBrec + self.rightBoundBrec
        gradA = math.fabs(self.rightBoundArec - self.leftBoundArec)
        gradB = math.fabs(self.rightBoundBrec - self.leftBoundBrec)
        if self.decisiontype == "naive":
            ratio = float(interA)/float(interB+interA)
            self.Brec = math.floor(self.MaxReceptors*ratio)
            self.Arec = self.MaxReceptors - self.Brec
        elif self.decisiontype == "measured":
            ratio1 = float(interA) / float(interB + interA)
            ratio2 = float(exterA)/float(exterB+exterA)
            self.Brec = math.floor(self.MaxReceptors * ((ratio1+ratio2)*(1.0/2.0)))
            self.Arec = self.MaxReceptors - self.Brec
        elif self.decisiontype == "prediction":
            ratio1 = float(interA) / float(interB + interA)
            ratio2 = float(exterA) / float(exterB + exterA)
            ratio3 = float(gradB) /float(gradA + gradB)
            self.Brec = math.floor(self.MaxReceptors * ((ratio1 + ratio2 + ratio3) * (1.0 / 3.0)))
            self.Arec = self.MaxReceptors - self.Brec
        elif self.decisiontype == "non":
            return
        return

    def divideStats(self):
        #mutation: either doesn't mutate or substracts r adds one receptor
        if self.mutate:
            newArec = math.floor(self.Arec) + self.newRand.getNewInt(0, 21) - 10
            if newArec < 0:
                newArec = 0
            newBrec = self.MaxReceptors - newArec
        else:
            newArec = self.Arec
            newBrec = self.Brec
        return [newArec, newBrec, self.MaxReceptors, math.floor(self.Amol/2), math.floor(self.Bmol/2), math.floor(self.ATP/2), math.floor(self.biomass/2), self.generation+1, self.distanceTrav]

    def toString(self):
        return str(self.Arec) +", "+ str(self.Brec) +", "+ str(self.Amol) +", "+ str(self.Bmol) + ", " + str(self.ATP) +", "+ str(self.biomass) + ", " + str(self.generation) +", "+ str(self.distanceTrav) + " ID: "+ str(self.ID)

    def toString2(self):
        return str(self.Arec) +", "+ str(self.Brec)