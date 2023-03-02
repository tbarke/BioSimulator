import math
import utils
import numpy as np
import log
l = log.log()

class cell(object):
#cell attributes: A receptors, B receptors, A molecules, B molecules, Bound A receptors, and Bound B receptors
    def __init__(self, newRand, Arec, Brec, Amol, Bmol, config):
        self.config = config
        self.Arec = Arec
        self.Brec = Brec
        self.MaxReceptors = config.cellStats.maxRec
        self.Amol = Amol
        self.Bmol = Bmol
        self.ATP = config.cellStats.ATP
        self.generation = config.cellStats.generation
        self.distanceTrav = config.cellStats.distTrav
        self.biomass = config.cellStats.biomass
        self.AbsorbtionRate = config.cellMetaStats.absorptionRate
        self.ReceptorConsumptionRate = config.cellMetaStats.receptorConsumptionRate
        self.survivalCost = config.cellMetaStats.survivalCost
        self.VelocityMultiplier = config.cellMetaStats.velocityMultiplier
        self.noise = config.cellMetaStats.noise
        self.receptor_mode = config.cellMetaStats.receptorMode
        self.Combined_portion = config.cellMetaStats.combinedPortion
        self.Divide_portion = config.cellMetaStats.dividePortion
        self.adaptive_ratio = config.cellMetaStats.adaptiveRatio
        self.mutate = True
        if config.cellMetaStats.mutate == "non":
            self.mutate = False
        self.decisiontype = config.cellMetaStats.decisiontype

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

        self.cell_stress_level = 0
        self.adaptive_ratio = 0

        self.dissocociation = config.cellMetaStats.dissocociationConstant


    def Absorb(self, AconLeft, AconRight, BconLeft, BconRight, timestep):
        retA = utils.calculateAbsorbtion((AconLeft + AconRight) / 2, self.AbsorbtionRate*timestep)
        retB = utils.calculateAbsorbtion((BconLeft + BconRight) / 2, self.AbsorbtionRate*timestep)
        #print(retA)
        self.Amol = self.Amol + retA
        self.Bmol = self.Bmol + retB
        #complex = utils.ATPcomplex(self.Amol, self.Bmol, 1, 1)
        #self.ATP += complex[0]
        #self.Amol -= complex[1]
        #self.Bmol -= complex[2]
        return [retA, retB]

    #returns the velocity of the cell
    def velocity(self, AconLeft, AconRight, BconLeft, BconRight, alphaA, alphaB, timeStep):
        leftArec = math.floor(self.Arec/2)
        rightArec = leftArec
        leftBrec = math.floor(self.Brec/2)
        rightBrec = leftBrec
        self.AconLeft = AconLeft
        self.AconRight = AconRight
        self.BconLeft = BconLeft
        self.BconRight = BconRight
        self.leftBoundArec = utils.calculateBoundKinetic(leftArec, AconLeft, self.dissocociation, self.newRand, self.noise, self.receptor_mode)
        self.rightBoundArec = utils.calculateBoundKinetic(rightArec, AconRight, self.dissocociation, self.newRand, self.noise, self.receptor_mode)
        self.leftBoundBrec = utils.calculateBoundKinetic(leftBrec, BconLeft, self.dissocociation, self.newRand, self.noise, self.receptor_mode)
        self.rightBoundBrec = utils.calculateBoundKinetic(rightBrec, BconRight, self.dissocociation, self.newRand, self.noise, self.receptor_mode)

        A_ratio = 0.5
        if self.Amol + self.Bmol > 0:
            A_ratio = self.Amol/(self.Amol + self.Bmol)

        adjusted_rightArec = self.rightBoundArec*(1-A_ratio)
        adjusted_rightBrec = self.rightBoundBrec*A_ratio
        adjusted_leftArec = self.leftBoundArec*(1-A_ratio)
        adjusted_leftBrec = self.leftBoundBrec*A_ratio

        wrong_rightArec = self.rightBoundArec * (A_ratio)
        wrong_rightBrec = self.rightBoundBrec * (1-A_ratio)
        wrong_leftArec = self.leftBoundArec * A_ratio
        wrong_leftBrec = self.leftBoundBrec * (1-A_ratio)
        if self.decisiontype == "adjusted":
            vel = 2*((adjusted_rightArec + adjusted_rightBrec) - (adjusted_leftArec+adjusted_leftBrec))
        if self.decisiontype == "drastic2":
            if self.Amol < self.Bmol:
                vel = 2*(self.rightBoundArec - self.leftBoundArec)
            else:
                vel = 2*(self.rightBoundBrec - self.leftBoundBrec)
        elif self.decisiontype == "wrong":
            vel = 2 *((wrong_rightArec + wrong_rightBrec) - (wrong_leftArec + wrong_leftBrec))
        elif self.decisiontype == "ratio":
            vel = 20 * ((self.rightBoundBrec - self.leftBoundBrec) * self.config.cellMetaStats.AratioInt * A_ratio) + ((self.rightBoundArec - self.leftBoundArec) * (1-self.config.cellMetaStats.AratioInt)* (1-A_ratio))
        #rightBound = utils.calculateBoundKinetic(rightBrec+rightArec, BconRight+AconRight, dissocociation, self.newRand, self.noise, self.receptor_mode)
        #leftBound = utils.calculateBoundKinetic(leftBrec+leftArec, BconLeft+AconLeft, dissocociation, self.newRand, self.noise, self.receptor_mode)

        #print(self.leftBoundArec)
        #print(self.leftBoundBrec)
        #print(self.rightBoundArec)
        #print(self.rightBoundBrec)

        #adopting right positive notation
        elif self.decisiontype == "non" or self.decisiontype == "measured" or self.decisiontype == "counter" or self.decisiontype == "drastic":
            vel = ((self.rightBoundArec + self.rightBoundBrec) - (self.leftBoundArec + self.leftBoundBrec))
        #vel = math.floor((rightBound -leftBound) * self.VelocityMultiplier)
        #vel = math.floor(((self.leftBoundArec + self.leftBoundBrec) - (self.rightBoundArec + self.rightBoundBrec)) * self.VelocityMultiplier)
        #print(vel)
        #exit(-1)
        self.vel = vel
        vel = vel*(1/self.config.simParams.locationStep)*timeStep
        rand = np.random.randint(-50, 50, 1)[0]
        #print(rand)
        #if self.vel > 10:
        #   self.vel = 10
        #elif self.vel < -10:
        #    self.vel = -10
        return vel #math.floor(((self.rightBoundArec + self.rightBoundBrec) - (self.leftBoundArec + self.leftBoundBrec))*self.VelocityMultiplier)

    def consumption(self, velocity, fullDie, timeStep):
        AllCost = self.survivalCost
        combined_ratio = self.Combined_portion

        #if self.Amol + self.Bmol != 0:

        #ratio = self.Amol/(self.Bmol+self.Amol)
        #self.Amol -= math.floor((2*AllCost*(1-combined_ratio)*ratio)*timeStep)
        #self.Bmol -= math.floor((2*AllCost*(1-combined_ratio)*(1-ratio))*timeStep)
        self.Amol -= math.floor((AllCost*combined_ratio)*timeStep)
        self.Bmol -= math.floor((AllCost*combined_ratio)*timeStep)

        if self.Amol < self.Bmol:
            self.cell_stress_level = self.Amol
            if self.Amol > self.survivalCost*5:
                self.cell_stress_level = self.survivalCost*5
        else:
            self.cell_stress_level = self.Bmol
            if self.Bmol > self.survivalCost*5:
                self.cell_stress_level = self.survivalCost*5

        self.cell_stress_level = (self.survivalCost*5 - self.cell_stress_level)/(self.survivalCost*5)


        if self.Amol < 0 or self.Bmol < 0:
            if fullDie:
                return False
        return True

        if self.Amol + self.Bmol <= 0:
            self.Amol = 0
            self.Bmol = 0
            if fullDie:
                return False
        if self.Amol < 0:
            self.Bmol += self.Amol
            self.Amol = 0
        elif self.Bmol < 0:
            self.Amol += self.Bmol
            self.Bmol = 0
        return True

    def incrementDis(self, distance):
        self.distanceTrav += distance

    def divide(self):
        self.biomass += 1
        if self.Amol > self.survivalCost*5 and self.Bmol > self.survivalCost*5:
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
        if self.decisiontype == "drastic":
            if interA > interB:
                self.Brec = self.MaxReceptors
                self.Arec = 0
            else:
                self.Arec = self.MaxReceptors
                self.Brec = 0
        if self.decisiontype == "measured":
            #ratio2 = 0.999
            if float(interA+interB) != 0:
                ratio = float(interA)/float(interB+interA)
            else:
                ratio = 0.5
            #self.Brec = math.floor(self.MaxReceptors*ratio*(1-self.adaptive_ratio))+math.floor(self.MaxReceptors*0.5*self.adaptive_ratio)
            self.Brec = math.floor(self.MaxReceptors*ratio)
            self.Arec = self.MaxReceptors - self.Brec
        if self.decisiontype == "ratio":
            self.Brec = math.floor(self.MaxReceptors*self.config.cellMetaStats.Aratio)
            self.Arec = self.MaxReceptors - self.Brec
        if self.decisiontype == "counter":
            if float(interA+interB) != 0:
                ratio = float(interA)/float(interB+interA)
            else:
                ratio = 0.5
            self.Arec = math.floor(self.MaxReceptors * ratio)
            self.Brec = self.MaxReceptors - self.Brec

        elif self.decisiontype == "measured2":
            ratio1 = float(interA) / float(interB + interA)
            ratio2 = float(exterA)/float(exterB+exterA)
            self.Brec = math.floor(self.MaxReceptors * ((ratio1+ratio2)*(1.0/2.0)))
            self.Arec = self.MaxReceptors - self.Brec
        elif self.decisiontype == "prediction":
            if float(interA+interB) != 0:
                ratio = float(interA)/float(interB+interA)
            else:
                ratio = 0
            #print(self.cell_stress_level)
            Bratio = self.Brec / (self.Brec + self.Arec)
            self.Brec = math.floor((1-self.cell_stress_level)*self.MaxReceptors*Bratio) + math.floor((self.cell_stress_level)*self.MaxReceptors*ratio)
            self.Arec = self.MaxReceptors - self.Brec

            #ratio1 = float(interA) / float(interB + interA)
            #ratio2 = float(exterA) / float(exterB + exterA)
            #ratio3 = float(gradB) /float(gradA + gradB)
            #self.Brec = math.floor(self.MaxReceptors * ((ratio1 + ratio2 + ratio3) * (1.0 / 3.0)))
            #self.Arec = self.MaxReceptors - self.Brec
        elif self.decisiontype == "non" or self.decisiontype == "adjusted" or self.decisiontype == "wrong" or self.decisiontype == "drastic2":
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