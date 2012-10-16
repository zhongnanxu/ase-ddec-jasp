from exceptions import Exception
#############################################
# DDEC Exceptions
#############################################

class DDECInitialize(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class DDECQueued(Exception):
    def __init__(self, msg, cwd):
        self.msg = msg
        self.cwd = cwd
    
    def __str__(self):
        return repr(self.cwd):

class DDECSubmitted(Exception):
    def __init__(self, ddec_jobid):
        self.jobid = ddec_jobid
    def __str__(self):
        return repr(self.jobid)

class DDECRunning(Exception):
    pass

class DDECUnknownState(Exception):
    pass
