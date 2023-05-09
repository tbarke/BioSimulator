import testBench
import log
import utils
import traceback
import io

dateTime = str(utils.getTodaysDate()) + '_' + utils.getTime()
l = log.log('Logs/'+dateTime + '.log')

try:
    path = testBench.createNewConfig('vonMisesLargeDissIdentical_05-05-23.cfg')
    l.log(path)
    testBench.main()
except Exception as e:
    l.handleException(e)

l.exit()
