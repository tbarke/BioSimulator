import io
import traceback

class log:
    #global variables
    logfile = ''
    currentMessage = ''

    def __init__(self, logFile = None):
        if logFile:
            log.logfile = logFile
            with open(logFile, 'w') as file:
                pass

    def saveMessage(self):
        with open(log.logfile, 'a') as file:
            file.write(log.currentMessage)
        log.currentMessage = ''

    def exit(self, code = None):
        self.saveMessage()
        if code:
            exit(code)
        else:
            exit()

    def log(self, message, prin = True):
        string_buffer = io.StringIO()
        print(message, file=string_buffer)
        if prin:
            print(message)
        log.currentMessage += string_buffer.getvalue()
        if len(log.currentMessage) > 1000:
            self.saveMessage()

    def handleException(self, e):
        try:
            raise e
        except Exception:
            string_buffer = io.StringIO()
            print(traceback.format_exc(), file=string_buffer)
            self.log(string_buffer.getvalue())

