from collections import OrderedDict

from .DefaultCard import DefaultCardInterface
from .dihiggs import DiHiggsInterface
from .eHDECAY import eHDECAYInterface
from .EWPO import EWPOInterface
from .Lilith import LilithInterface
from .SignalStrengths import SignalStrengthsInterface
from .Translate import TranslateInterface

__all__ = [
    "DefaultCardInterface",
    "DiHiggsInterface",
    "eHDECAYInterface",
    "EWPOInterface",
    "LilithInterface",
    "SignalStrengthsInterface",
    "TranslateInterface"
]

_all_interfaces = OrderedDict({k: globals()[k] for k in __all__})
