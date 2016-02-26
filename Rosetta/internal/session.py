import sys
import settings
from StringIO import StringIO

__once = set()
output = StringIO()
suppressed = StringIO()

def stdout(msg, log=False):
    '''Display a message.'''
    print >> output, msg
    print msg
    
def log(msg):
    '''Display a message if settings.silent == False.'''
    
    if not settings.silent: 
        stdout(msg)
    else:
        print >> suppressed, msg
    

def once(msg):
    '''
    Display a message that should be shown only once during runtime. 
    Once printed, msg is stored in module state variable '__once'.
    '''
    if msg not in __once:
        log(msg)
        __once.add(msg)

def verbose(msg):
    '''
    Display a message that should only be shown if settings.verbose == True 
    and settings.silent == False.
    '''
    if settings.verbose: 
        log(msg)
    else:
        print >> suppressed, msg
    
def query(question, default="yes"):
    """
    Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    Recipe borrowed from http://code.activestate.com/recipes/577058/
    """
    valid = {"yes":True, "y":True, "ye":True,
             "no":False, "n":False}
    
    if settings.silent: 
        print >> suppressed, question+' [{}]'.format(default)
        return valid[default]
    elif settings.force: 
        stdout(question+' [{}]'.format(default))
        return valid[default]
    
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        log(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            log("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

def exit(code=0):
    '''
    Exit rosetta session by raisin sys.exit(code). The contents of 
    settings.output and settings.suppressed are written to 'rosetta.log' and 
    'rosetta.suppressed.log' respectively.
    '''
    with open('rosetta.log', 'w') as logfile:
        logfile.write(output.getvalue())
    with open('rosetta.suppressed.log', 'w') as logfile:
        logfile.write(suppressed.getvalue())
    log('Exit.')
    log('Output and suppressed output written to rosetta.log and '
           'rosetta.suppressed.log respectively.')
    sys.exit(code)
    