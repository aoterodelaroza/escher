#!/usr/bin/env python

import time
import cProfile

def profile_this(fn):
    '''
    A decorator for profiling functions
    @profile_this
    def foo():
    '''

    def profiled_fn(*args, **kwargs):
        '''
        .profile can be seen with RunSnakeRun
        '''
        fpath = fn.__name__ + '.profile'
        prof = cProfile.Profile()
        ret = prof.runcall(fn, *args, **kwargs)
        prof.dump_stats(fpath)
        return ret

    return profiled_fn


def time_this(fn):
    def timed_fn(*args, **kwargs):
        start_s = time.time()
        ret = fn(*args, **kwargs)
        elapsed_s = time.time() - start_s
        print '{} took {} s'.format(fn.__name__, elapsed_s)
        return ret
    return timed_fn
        
