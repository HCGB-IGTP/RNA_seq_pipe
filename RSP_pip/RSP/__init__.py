from pkg_resources import get_distribution
try:
    __version__ = get_distribution('RSP').version
except:
    try:
        __version__ = pkg_resources.resorce_filename('RSP', VERSION)
    except:
        __version__ = 'local'

__all__ = [
    'modules',
    'data',
    'scripts',
    'config'
    ]

from RSP import *
