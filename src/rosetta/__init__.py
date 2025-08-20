__author__ = ('Adam Falkowski, Benjamin Fuks, Kentarou Mawatari, Ken Mimasu, '
              'Francesco Riva & Veronica Sanz')
__date__ = '(undefined)'
__url__ = 'http://rosetta.hepforge.org'


from .internal.machinery import bases as implemented_bases

from .internal import session, settings

from .bases import *
from ._version import __version__