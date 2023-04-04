import configparser
import math
import re
import log
l = log.log()

class configuration(object):

    class ConcParams:
        def __init__(self, simParams):
            self.concProfile = "vonMises"
            self.constantConc = 10
            self.uphillHigh = 14
            self.uphillLow = 0
            self.diffCoeff = 5
            self.arp = "periodic"
            self.repeatFrequency = 0
            self.magnitude = 100000.0
            self.locationA = 25
            self.locationB = 75
            self.locationC = 50
            self.concDiffuse = "static"
            self.VonMisesVar = 10 * ((2*math.pi)/simParams.length)
            self.VonMisesMagnitude = 175
            self.VonMisesAOffset = -25 * (1/simParams.locationStep)
            self.VonMisesBOffset = 25 * (1/simParams.locationStep)
            self.VonMisesCOffset = 0 * (1 / simParams.locationStep)

    class SimParams:
        def __init__(self):
            self.length = 100
            self.locationStep = 0.1
            self.simTimeStep = 0.01
            self.simLength = 20
            self.presetA_time = []
            self.presetB_time = []
            self.presetA_loc = []
            self.presetB_loc = []
            self.AconcFile = 'None'
            self.BconcFile = 'None'
            self.presetSave = False
            self.presetBool = False
            self.uniformConcStart = False
            self.simulationConcStart = True
            self.lam = 0.4
            self.gamma = 0.1
            self.aknoght = 200
            self.highSpace = 2000
            self.lowSpace = 1000

    class CellMetaStats:
        def __init__(self):
            self.stress = 0.01
            self.absorptionRate = 1/self.stress
            self.receptorConsumptionRate = 0.0
            self.survivalCost = 250
            self.velocityMultiplier = 1.0
            self.mutate = "non"
            self.decisiontype = "non"
            self.noise = 0.0
            self.receptorMode = "simulate"
            self.combinedPortion = 1
            self.dividePortion = 1
            self.adaptiveRatio = 0
            self.cellLocations = []
            self.fullDivide = True
            self.fullDie = True
            self.dissocociationConstant = 2.0
            self.Aratio = 0.5
            self.AratioInt = 0.5

    class OutputFlags:
        def __init__(self):
            self.totalcells = True
            self.totalReceptors = True
            self.boundReceptors = True
            self.internalAB = True
            self.totalConcs = True
            self.cellLocations = True
            self.cellMovement = True
            self.concProfile = True

    class RunOutputFlags:
        def __init__(self):
            self.save = True
            self.compressSave = False
            self.totalMI = True
            self.averageMI = True
            self.compositeMI = True
            self.weightedMI = True
            self.averageGrowth = True
            self.SIAsVar = True
            self.SIAsEntropy = True

    class CellStats:
        def __init__(self):
            self.Arec = 200
            self.Brec = 200
            self.maxRec = 400
            self.Amol = 0
            self.Bmol = 0
            self.ATP = 0
            self.biomass = 5
            self.generation = 1
            self.distTrav = 0
            self.startID = 0

    class RunStats:
        def __init__(self):
            self.stressArray = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            self.noiseArray = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            self.enviornmentArray = ["vonMises", "constUpDown", "constant", "manaFromHeaven"]
            self.cellStrategiesArray = ["non", "measured", "counter", "adjusted", "drastic", "drastic2"]
            self.cellRatioAEmphasis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            self.cellRatioAIntEmphasis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            self.cellRatioAEmphasisFlag = True
            self.cellRatioAIntEmphasisFlag = True
            self.cellStrategiesArrayFlag = True
            self.enviornmentArrayFlag = True
            self.runStress = True
            self.runNoise = False
            self.runDetermine = False
            self.beginningRandInt = 1
            self.save = True
            self.saveDir = "Data"

    def __init__(self):
        self.simParams = self.SimParams()
        self.concParams = self.ConcParams(self.simParams)
        self.cellMetaStats = self.CellMetaStats()
        self.outputFlags = self.OutputFlags()
        self.runOutputFlags = self.RunOutputFlags()
        self.runStats = self.RunStats()
        self.cellStats = self.CellStats()

    def parse_string(self, input_str):
        def parse_value(input_string):
            # Try to parse the string as an integer
            match = re.match(r'^[-+]?\d+$', input_string)
            if match:
                return int(input_string)

            # Try to parse the string as a float
            match = re.match(r'^[-+]?\d+\.\d+$', input_string)
            if match:
                return float(input_string)

            # Try to parse the string as a boolean
            if input_string.lower() in ['true', 'false']:
                return input_string == 'True' or input_string == 'true'


            # Return the string as is
            if input_string[0] == '\'' and input_string[len(input_string)-1] == '\'':
                input_string = input_string.replace('\'','')

            return input_string

        def parse_string_list(input_string):
            # Check if the string starts and ends with square brackets
            if not input_string.startswith('[') or not input_string.endswith(']'):
                return None

            # Remove square brackets from the beginning and end of the string
            input_string = input_string[1:-1]

            if input_string.strip() == "":
                return []

            # Split the string by commas and strip whitespace from each element
            string_list = [elem.strip() for elem in input_string.split(',')]

            # Check if any element in the list is empty or contains square brackets
            for elem in string_list:
                if not elem or '[' in elem or ']' in elem:
                    return None

            # Return the list of strings
            return string_list

        input_str = input_str.strip()
        list = parse_string_list(input_str)
        if list == []:
            return list
        if list:
            ret_list = []
            for l in list:
                ret_list.append(parse_value(l))
            return ret_list

        return parse_value(input_str)

    """""
    def parse_string(self, input_str):
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
    """""
    def readConfig(self, filename, supressWarnings = True):
        config = configparser.ConfigParser()
        config.read(filename)

        mapNames = {'cellMetaStats': 'CellMetaStats', 'cellStats': 'CellStats', 'concParams': 'ConcParams', 'outputFlags': 'OutputFlags', 'runOutputFlags': 'RunOutputFlags', 'runStats': 'RunStats','simParams': 'SimParams'}
        paramOrg = self.__dict__
        for key1, value1 in paramOrg.items():
            paramnames = value1.__dict__
            for key2 in paramnames:
                try:
                    paramnames[key2] = self.parse_string(config.get(mapNames[key1], key2))
                except configparser.NoOptionError:
                    if not supressWarnings:
                        l.log("Warning: could not load in: " + mapNames[key1] + "." + key2 + ", Probably loading in from an old version config. Using default Value.")

    def writeConfig(self, filename):
        config = configparser.ConfigParser()
        config['RunStats'] = {'stressArray': self.runStats.stressArray,
                              'noiseArray': self.runStats.noiseArray,
                              'enviornmentArray': self.runStats.enviornmentArray,
                              'cellStrategiesArray': self.runStats.cellStrategiesArray,
                              'runStress': self.runStats.runStress,
                              'cellStrategiesArrayFlag': self.runStats.cellStrategiesArrayFlag,
                              'enviornmentArrayFlag': self.runStats.enviornmentArrayFlag,
                              'runNoise': self.runStats.runNoise,
                              'runDetermine': self.runStats.runDetermine,
                              'beginningRandInt': self.runStats.beginningRandInt,
                              'cellRatioAEmphasis': self.runStats.cellRatioAEmphasis,
                              'cellRatioAEmphasisFlag': self.runStats.cellRatioAEmphasisFlag,
                              'cellRatioAIntEmphasis': self.runStats.cellRatioAIntEmphasis,
                              'cellRatioAIntEmphasisFlag': self.runStats.cellRatioAIntEmphasisFlag,
                              'save': self.runStats.save,
                              'saveDir' : self.runStats.saveDir
        }

        config['CellStats'] = {'Arec': self.cellStats.Arec,
                              'Brec': self.cellStats.Brec,
                              'maxRec': self.cellStats.maxRec,
                              'Amol': self.cellStats.Amol,
                              'Bmol': self.cellStats.Bmol,
                              'ATP': self.cellStats.ATP,
                              'biomass': self.cellStats.biomass,
                                'generation': self.cellStats.generation,
                               'distTrav': self.cellStats.distTrav,
                               'startID': self.cellStats.startID
                              }

        config['RunOutputFlags'] = {'save': self.runOutputFlags.save,
                               'compressSave': self.runOutputFlags.compressSave,
                               'totalMI': self.runOutputFlags.totalMI,
                               'averageMI': self.runOutputFlags.averageMI,
                               'compositeMI': self.runOutputFlags.compositeMI,
                               'weightedMI': self.runOutputFlags.weightedMI,
                               'averageGrowth': self.runOutputFlags.averageGrowth,
                               'SIAsVar': self.runOutputFlags.SIAsVar,
                               'SIAsEntropy': self.runOutputFlags.SIAsEntropy
                               }

        config['OutputFlags'] = {'totalcells': self.outputFlags.totalcells,
                                    'totalReceptors': self.outputFlags.totalReceptors,
                                    'boundReceptors': self.outputFlags.boundReceptors,
                                    'internalAB': self.outputFlags.internalAB,
                                    'totalConcs': self.outputFlags.totalConcs,
                                    'cellLocations': self.outputFlags.cellLocations,
                                    'cellMovement': self.outputFlags.cellMovement,
                                    'concProfile': self.outputFlags.concProfile
                                    }

        config['CellMetaStats'] = {'stress': self.cellMetaStats.stress,
                                    'absorptionRate': self.cellMetaStats.absorptionRate,
                                    'receptorConsumptionRate': self.cellMetaStats.receptorConsumptionRate,
                                    'survivalCost': self.cellMetaStats.survivalCost,
                                    'velocityMultiplier': self.cellMetaStats.velocityMultiplier,
                                    'mutate': self.cellMetaStats.mutate,
                                    'decisiontype': self.cellMetaStats.decisiontype,
                                    'noise': self.cellMetaStats.noise,
                                    'receptorMode': self.cellMetaStats.receptorMode,
                                    'combinedPortion': self.cellMetaStats.combinedPortion,
                                    'dividePortion': self.cellMetaStats.dividePortion,
                                    'adaptiveRatio': self.cellMetaStats.adaptiveRatio,
                                    'cellLocations': self.cellMetaStats.cellLocations,
                                    'fullDivide': self.cellMetaStats.fullDivide,
                                    'fullDie': self.cellMetaStats.fullDie,
                                    'Aratio': self.cellMetaStats.Aratio,
                                    'AratioInt': self.cellMetaStats.AratioInt,
                                    'dissocociationConstant': self.cellMetaStats.dissocociationConstant
                                    }

        config['SimParams'] = {'length': self.simParams.length,
                                   'locationStep': self.simParams.locationStep,
                                   'simTimeStep': self.simParams.simTimeStep,
                                   'simLength': self.simParams.simLength,
                                   'presetA_time': self.simParams.presetA_time,
                                   'presetB_time': self.simParams.presetB_time,
                                   'presetA_loc': self.simParams.presetA_loc,
                                   'presetB_loc': self.simParams.presetB_loc,
                                   'presetBool': self.simParams.presetBool,
                                   'uniformConcStart': self.simParams.uniformConcStart,
                                   'simulationConcStart': self.simParams.simulationConcStart,
                                   'lam': self.simParams.lam,
                                   'gamma': self.simParams.gamma,
                                   'aknoght': self.simParams.aknoght,
                                   'presetSave': self.simParams.presetSave,
                                   'highSpace': self.simParams.highSpace,
                                   'lowSpace': self.simParams.lowSpace,
                                   'AconcFile': self.simParams.AconcFile,
                                   'BconcFile': self.simParams.BconcFile

                               }

        config['ConcParams'] = {'diffCoeff': self.concParams.diffCoeff,
                               'arp': self.concParams.arp,
                               'concProfile': self.concParams.concProfile,
                               'constantConc': self.concParams.constantConc,
                               'uphillHigh': self.concParams.uphillHigh,
                               'uphillLow': self.concParams.uphillLow,
                               'repeatFrequency': self.concParams.repeatFrequency,
                               'magnitude': self.concParams.magnitude,
                               'locationA': self.concParams.locationA,
                               'locationB': self.concParams.locationB,
                               'locationC': self.concParams.locationC,
                               'concDiffuse': self.concParams.concDiffuse,
                               'VonMisesVar': self.concParams.VonMisesVar,
                               'VonMisesMagnitude': self.concParams.VonMisesMagnitude,
                               'VonMisesAOffset': self.concParams.VonMisesAOffset,
                               'VonMisesBOffset': self.concParams.VonMisesBOffset,
                               'VonMisesCOffset': self.concParams.VonMisesCOffset
                               }
        with open(filename, 'w+') as configfile:
            config.write(configfile)
        return filename
