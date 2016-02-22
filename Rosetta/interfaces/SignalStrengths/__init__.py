from errors import SignalStrengthsImportError, SignalStrengthsInterfaceError
from ..eHDECAY.errors import eHDECAYImportError

try:
    from ..eHDECAY.eHDECAY import SM_BR
    use_eHDECAY = True
except eHDECAYImportError:
    use_eHDECAY = False
    