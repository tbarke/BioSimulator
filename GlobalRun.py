import testBench
import log
import utils
import traceback
import io

dateTime = str(utils.getTodaysDate()) + '_' + utils.getTime()
l = log.log('Logs/'+dateTime + '.log')

try:
    testBench.main()
except Exception as e:
    l.handleException(e)

l.exit()
