from ..errors import TranslationError

__doc__ = '''
The translation decorator used to tag a translation function in a Rosetta basis 
implementation.
'''
################################################################################
def parametrised(dec):
    '''
    Decorator for a decorator allowing it to take arguments
    '''
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)
        return repl
    return layer

# translation decorator to identify translation functions and their target as 
# well as handle errors.
@parametrised
def translation(func, target):
    '''
    translation decorator for Rosetta to identify translation functions and 
    their target basis
    '''
    def labelled(*args, **kwargs):
        try:
            return func(*args,**kwargs)
        except Exception as e:
            import traceback
            raise TranslationError(traceback.format_exc())

    labelled.__name__ = func.__name__
    labelled._target = target

    return labelled


@parametrised
def translation_from(func, source):
    '''
    translation decorator for Rosetta to identify translation function from an 
    external basis implementation, 'source'
    '''
    def labelled(*args, **kwargs):
        try:
            return func(*args,**kwargs)
        except Exception as e:
            import traceback
            raise TranslationError(traceback.format_exc())

    labelled.__name__ = func.__name__
    labelled._source = source

    return labelled

def derived_input(func):
    '''
    translation decorator for Rosetta to identify functions that derive input 
    parameters from others
    '''
    def labelled(*args, **kwargs):
        try:
            return func(*args,**kwargs)
        except Exception as e:
            import traceback
            raise TranslationError(traceback.format_exc())

    labelled.__name__ = func.__name__
    labelled._derived_input = func.__name__

    return labelled

################################################################################
