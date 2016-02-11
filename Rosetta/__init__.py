import os

files = os.listdir(os.path.dirname(__file__))

__all__ = [f.replace('.py','') for f in files
         if '.py'==f[-3:] and f!='__init__.py']

