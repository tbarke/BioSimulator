import configparser
import math

class output(object):

    class PrimativeOutput:
        def __init__(self, totalcellsFile = None, totalReceptorsFile = None, boundReceptorsFile = None, internalABFile = None, totalConcsFile = None, cellLocationsFile = None, cellMovementFile = None, enviornmentConcsFile = None):
            self.totalcells = totalcellsFile
            self.totalReceptors = totalReceptorsFile
            self.boundReceptors = boundReceptorsFile
            self.internalAB = internalABFile
            self.totalConcs = totalConcsFile
            self.cellLocationsFile = cellLocationsFile
            self.cellMovementFile = cellMovementFile
            self.enviornmentConcsFile = enviornmentConcsFile

    class CalcOutput:
        def __init__(self, totalMIFile = None, averageMIFile = None, compositeMIFile = None, weightedMIFile = None, averageGrowthFile = None, SIAsVarFile = None, SIAsEntropyFile = None):
            self.totalMI = totalMIFile
            self.averageMI = averageMIFile
            self.compositeMI = compositeMIFile
            self.weightedMI = weightedMIFile
            self.averageGrowth = averageGrowthFile
            self.SIAsVar = SIAsVarFile
            self.SIAsEntropy = SIAsEntropyFile

    def __init__(self, runName = None, totalcellsFile = None, totalReceptorsFile = None, boundReceptorsFile = None, internalABFile = None, totalConcsFile = None, cellLocationsFile = None, cellMovementFile = None, enviornmentConcsFile = None, totalMIFile = None, averageMIFile = None, compositeMIFile = None, weightedMIFile = None, averageGrowthFile = None, SIAsVarFile = None, SIAsEntropyFile = None):
        self.runName = runName
        self.primativeOutput = self.PrimativeOutput(totalcellsFile,totalReceptorsFile, boundReceptorsFile, internalABFile, totalConcsFile, cellLocationsFile, cellMovementFile, enviornmentConcsFile)
        self.calcOutput = self.CalcOutput(totalMIFile, averageMIFile, compositeMIFile, weightedMIFile, averageGrowthFile, SIAsVarFile, SIAsEntropyFile)