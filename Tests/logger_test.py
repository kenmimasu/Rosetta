import sys

sys.path.append('../')

from logger import Logger

class A(object):
    def __init__(self):
        self.log = Logger()
    
monkeys = A()
bears = B()

monkeys.log.onetime('malarkey')
bears.log.onetime('malarkey')