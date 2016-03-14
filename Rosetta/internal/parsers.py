from ..interfaces import *
from ..interfaces import __all__ as interface_names

interfaces = {n:globals()[n] for n in interface_names}

def set_subparsers(subparsers):
    parsers = {}
    for name, intr in interfaces.iteritems():
            
        parsers[name] = subparsers.add_parser(
                               intr.interface, 
                               help = intr.helpstr, 
                               description = intr.description
                           )
                                             
        for args, kwargs in intr.parser_args.iteritems():
            parsers[name].add_argument(*args, **kwargs)
    
        parsers[name].set_defaults(interface = intr)