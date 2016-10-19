# from ..interfaces import *
# from ..interfaces import __all__ as interface_names
from ..interfaces import _all_interfaces as interfaces
# print interfaces
# interfaces = {n:globals()[n] for n in interface_names}

# def set_subparsers(subparsers):
#     parsers = {}
#     for name, intr in interfaces.iteritems():
#
#         parsers[name] = subparsers.add_parser(
#                                intr.interface,
#                                help = intr.helpstr,
#                                description = intr.description
#                                )
#
#         for args, kwargs in intr.parser_args.iteritems():
#             parsers[name].add_argument(*args, **kwargs)
#
#         parsers[name].set_defaults(interface = intr)

def set_subparsers(subparsers):
    parsers = {}
    
    # Put translate interface first
    trans_item = interfaces.pop('TranslateInterface')
    interfaces_sorted = [('TranslateInterface', trans_item)]+interfaces.items()

    for name, intr in interfaces_sorted:
            
        parsers[name] = subparsers.add_parser(
                               intr.interface, 
                               help = intr.helpstr, 
                               description = intr.description
                               )
                                             
        for args, kwargs in intr.parser_args.iteritems():
            parsers[name].add_argument(*args, **kwargs)
    
        parsers[name].set_defaults(interface = intr)