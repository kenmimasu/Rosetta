from errors import SignalStrengthsImportError, SignalStrengthsInterfaceError

try:
    from ..eHDECAY.eHDECAY import SM_BR
    use_eHDECAY = True
except ImportError as e:
    use_eHDECAY = False
    
    
from interface import SignalStrengthsInterface
