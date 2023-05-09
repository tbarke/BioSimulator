import configparser
import math
import typing
import outputAnalysis
import log
l = log.log()

class output(object):
    class PrimitiveOutput:
        def __init__(self, totalcellsFile = None, totalReceptorsFile = None, boundReceptorsFile = None, internalABFile = None, totalConcsFile = None, cellLocationsFile = None, cellMovementFile = None, enviornmentConcsFile = None):
            self.totalcellsFile = totalcellsFile
            self.totalReceptorsFile = totalReceptorsFile
            self.boundReceptorsFile = boundReceptorsFile
            self.internalABFile = internalABFile
            self.totalConcsFile = totalConcsFile
            self.cellLocationsFile = cellLocationsFile
            self.cellMovementFile = cellMovementFile
            self.enviornmentConcsFile = enviornmentConcsFile

    class CalcOutput:
        def __init__(self, totalMIFile = None, averageMIFile = None, compositeMIFile = None, weightedMIFile = None, averageGrowthFile = None, SIAsVarFile = None, SIAsEntropyFile = None):
            self.totalMIFile = totalMIFile
            self.averageMIFile = averageMIFile
            self.compositeMIFile = compositeMIFile
            self.weightedMIFile = weightedMIFile
            self.averageGrowthFile = averageGrowthFile
            self.SIAsVarFile = SIAsVarFile
            self.SIAsEntropyFile = SIAsEntropyFile

    class RunClacOut:
        def __init__(self, MItrad = None, MI2D2D = None, MI2D1D = None, MI2D1Din = None, growth = None, intweightMImove = None,intABWeightind = None ):
            self.MItrad = MItrad
            self.MI2D2D = MI2D2D
            self.MI2D1D = MI2D1D
            self.MI2D1Din = MI2D1Din
            self.growth = growth
            self.intweightMImove = intweightMImove
            self.intABWeightind = intABWeightind

    def __init__(self, configFile = None, runName = None, totalcellsFile = None, totalReceptorsFile = None, boundReceptorsFile = None, internalABFile = None, totalConcsFile = None, cellLocationsFile = None, cellMovementFile = None, enviornmentConcsFile = None, totalMIFile = None, averageMIFile = None, compositeMIFile = None, weightedMIFile = None, averageGrowthFile = None, SIAsVarFile = None, SIAsEntropyFile = None):
        if runName:
            self.runName = runName
        else:
            self.runName = ''
        self.configFile = configFile
        self.fileName = ''
        self.PrimitiveOutput = self.PrimitiveOutput(totalcellsFile,totalReceptorsFile, boundReceptorsFile, internalABFile, totalConcsFile, cellLocationsFile, cellMovementFile, enviornmentConcsFile)
        self.CalcOutput = self.CalcOutput(totalMIFile, averageMIFile, compositeMIFile, weightedMIFile, averageGrowthFile, SIAsVarFile, SIAsEntropyFile)
        self.RunClacOut = self.RunClacOut()

    def refactorfile(self, file1, path):
        return path + '/' + file1

    def calculateMeasures(self, config, bins, path, MItradFlag = True, MI2D2DFlag = True, MI2D1DFlag = True, MI2D1DinFlag = True, growthFlag = True, intweightMImoveFlag = True, ABintdynamicWeightFlag = True):
        try:
            internalABFileRef = self.refactorfile(self.PrimitiveOutput.internalABFile, path)
            totalConcsFileRef = self.refactorfile(self.PrimitiveOutput.totalConcsFile, path)
            cellMovementFileRef = self.refactorfile(self.PrimitiveOutput.cellMovementFile, path)
            boundReceptorsFileRef = self.refactorfile(self.PrimitiveOutput.boundReceptorsFile, path)
            totalcellsFileRef = self.refactorfile(self.PrimitiveOutput.totalcellsFile, path)
            recsFileRef = self.refactorfile(self.PrimitiveOutput.totalReceptorsFile, path)

            MITrad, MI2D2D, MImove, growth, intweightMImove, intABWeightind = outputAnalysis.CalcData(config, bins, MITradFlag = MItradFlag, MI2d2dFlag = MI2D2DFlag, MIMoveFlag = MI2D1DFlag, growthFlag = growthFlag, intWeightFlag = intweightMImoveFlag, ABintdynamicWeightFlag = ABintdynamicWeightFlag, intABFile = internalABFileRef, extABFile = totalConcsFileRef, moveFile = cellMovementFileRef, boundFile = boundReceptorsFileRef, totalCellsFile =totalcellsFileRef, recsFile=recsFileRef)
            #TODO fix below and in write (object should be lower case to use the constructor, but this will break below)
            if MItradFlag:
                self.RunClacOut.MItrad = MITrad
            if MI2D2DFlag:
                self.RunClacOut.MI2D2D = MI2D2D
            if MI2D1DFlag:
                self.RunClacOut.MI2D1D = MImove
            if MI2D1DinFlag:
                pass
                #TODO implement
                #self.RunClacOut.MI2D1Din = MI2D1Din
            if growthFlag:
                self.RunClacOut.growth = growth
            if intweightMImoveFlag:
                self.RunClacOut.intweightMImove = intweightMImove
            if ABintdynamicWeightFlag:
                self.RunClacOut.intABWeightind = intABWeightind
        except Exception as e:
            l.handleException(e)

    def read(self, filename, supressWarnings = True):

        def parse_string(input_str):
            import re
            input_str = input_str.strip()
            # Check for integer
            if re.match(r'^[-+]?\d+$', input_str):
                return int(input_str)

            # Check for float
            if re.match(r'^[-+]?(\d+(\.\d*)?|\.\d+)$', input_str):
                return float(input_str)

            # Check for list of integers or floats
            if re.match(r'^\[(?:\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)\s*,?)*\]$', input_str):
                values_str = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)', input_str)
                values = [int(x) if x.isdigit() else float(x) for x in values_str]
                return values

            # Check for boolean
            if input_str.lower() == 'true':
                return True
            if input_str.lower() == 'false':
                return False

            # Default to returning the input string as a string
            return input_str

        config = configparser.ConfigParser()
        config.read(filename)

        paramOrg = self.__dict__
        for key1, value1 in paramOrg.items():
            if key1 == 'configFile' or key1 == 'runName' or key1 == 'fileName':
                paramOrg[key1] = config.get('nosection', key1)
                continue
            paramnames = value1.__dict__
            for key2 in paramnames:
                try:
                    paramnames[key2] = parse_string(config.get(key1, key2))
                except (configparser.NoOptionError, configparser.NoSectionError):
                    if not supressWarnings:
                        l.log("Warning: could not load in: " + key1 + "." + key2 + ", Probably loading in from an old version config. Using default Value.")

    def write(self, path, absolute = False):
        config = configparser.ConfigParser()

        priDict = {}
        calcDict = {}
        runClacDict = {}
        nosectionDict = {}
        paramOrg = self.__dict__
        primitive = (int, str, bool, float)
        for key1, value1 in paramOrg.items():
            if isinstance(value1, primitive):
                nosectionDict[key1] = value1
                continue
            paramnames = value1.__dict__
            for key2, value2 in paramnames.items():
                newValue2 = value2
                if not value2:
                    newValue2 = ''
                if key1 == "PrimitiveOutput":
                    priDict[key2] = newValue2
                elif key1 == "CalcOutput":
                    calcDict[key2] = newValue2
                elif key1 == "RunClacOut":
                    runClacDict[key2] = newValue2

        config['PrimitiveOutput'] = priDict
        config['CalcOutput'] = calcDict
        config['RunClacOut'] = runClacDict
        config['nosection'] = nosectionDict
        fullpath = path
        if not absolute:
            nameDir = self.runName.split('/')[0]
            fullpath = path + "/" + nameDir  + "/OutputProfiles/" + self.runName.replace("/","-") + ".cfg"
        with open(fullpath, 'w+') as configfile:
            config.write(configfile)

        return fullpath

