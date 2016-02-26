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

    # labelled.__name__ = func.__name__
    labelled._target = target

    return labelled
################################################################################
