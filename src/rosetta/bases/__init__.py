from rosetta.internal.basis.basis import Basis
from . import BSMCharacterisation
from . import HiggsBasis
from . import HiggsCharacterisation
from . import HiggsPO
from . import HISZ
from . import SILHBasis
from . import TemplateBasis
from . import WarsawBasis

__all__ = [
    "BSMCharacterisation",
    "HiggsBasis",
    "HiggsCharacterisation",
    "HiggsPO",
    "HISZ",
    "SILHBasis",
    "TemplateBasis",
    "WarsawBasis",
]

_all_bases: dict[str, type[Basis]] = {
    cls.name: cls
    for cls in map(lambda mod: getattr(globals()[mod], mod), __all__)
    if issubclass(cls, Basis)
}
