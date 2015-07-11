################################################################################
# decorator for a decorator to allow it to take an argument
def parametrised(dec):
    def layer(*args, **kwargs):
        def repl(f):
            return dec(f, *args, **kwargs)
        return repl
    return layer

# translation decorator to identify translation functions and their target
@parametrised
def translation(func, target):
    def labelled(*args,**kwargs):
        # print func.__dict__
        return func(*args,**kwargs)
    labelled._target = target
    return labelled
################################################################################
