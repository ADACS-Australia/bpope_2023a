# -*- coding: utf-8 -*-
from .. import config
from ..compat import Node, change_flags, theano, tt, is_tensor
import numpy as np
from functools import wraps
import logging
import sys

import time

import theano.d3viz as d3v
from theano.printing import pydotprint

logger = logging.getLogger("starry.ops")

__all__ = ["logger", "autocompile", "is_tensor", "clear_cache"]


booleans = (np.array(True).dtype,)
integers = (int, np.int16, np.int32, np.int64)
floats = (float, np.float16, np.float32, np.float64)


class CompileLogMessage:
    """
    Log a brief message saying what method is currently
    being compiled and print `Done` when finished.

    """

    def __init__(self, name, custom_message=None):
        self.name = name
        self.custom_message = custom_message
        self.locking = False

    def __enter__(self):
        if not config.message_lock:
            config.message_lock = True
            self.locking = True
            config.rootHandler.terminator = ""
            if self.custom_message is None:
                logger.info("Compiling `{0}`...".format(self.name))
            else:
                logger.info(self.custom_message)

    def __exit__(self, type, value, traceback):
        if self.locking:
            config.rootHandler.terminator = "\n"
            logger.info(" Done.")
            config.message_lock = False
            self.locking = False


def _get_type(arg):
    """
    Get the theano tensor type corresponding to `arg`.

    Note that arg must be one of the following:
        - a theano tensor
        - an integer (`int`, `np.int`, `np.int16`, `np.int32`, `np.int64`)
        - a numpy boolean (`np.array(True)`, `np.array(False)`)
        - a numpy float array with ndim equal to 0, 1, 2, or 3

    TODO: We could just do `tt.as_tensor_variable(arg).type` and then upcast...

    """
    ttype = type(arg)
    if is_tensor(arg):

        # Trivial
        return ttype

    else:

        # Cast lists to arrays
        if ttype in (list, tuple):
            if is_tensor(*arg):
                return type(arg[0])
            else:
                arg = np.array(arg)
                ttype = type(arg)

        # Determine the type
        if ttype in booleans:
            return tt.bscalar
        if ttype in integers:
            return tt.lscalar
        elif ttype in floats:
            return tt.dscalar
        elif hasattr(arg, "ndim"):
            if arg.ndim == 0:
                if arg.dtype in booleans:
                    return tt.bscalar
                elif arg.dtype in integers:
                    return tt.lscalar
                else:
                    return tt.dscalar
            elif arg.ndim == 1:
                if arg.dtype in booleans:
                    return tt.bvector
                elif arg.dtype in integers:
                    return tt.lvector
                else:
                    return tt.dvector
            elif arg.ndim == 2:
                if arg.dtype in booleans:
                    return tt.bmatrix
                elif arg.dtype in integers:
                    return tt.lmatrix
                else:
                    return tt.dmatrix
            elif arg.ndim == 3:
                if arg.dtype in booleans:
                    return tt.btensor3
                elif arg.dtype in integers:
                    return tt.ltensor3
                else:
                    return tt.dtensor3
            else:
                raise NotImplementedError(
                    "Invalid array dimension passed to @autocompile: {}.".format(
                        arg.ndim
                    )
                )
        else:
            raise NotImplementedError(
                "Invalid argument type passed to @autocompile: {}.".format(
                    ttype
                )
            )


def autocompile(func):
    """
    Wrap the method `func` and return a compiled version
    if none of the arguments are tensors.

    """

    @wraps(func)  # inherit docstring
    def wrapper(instance, *args):

        if is_tensor(*args):
            # print("is tensor: ", func.__name__, instance)
            # Just return the function as is
            # print("type: ", func.__name__, type(func(instance, *args)))
            return func(instance, *args)

        else:
            time.sleep(2)
            print("sleep: ", func.__name__)
            # Determine the argument types
            arg_types = tuple([_get_type(arg) for arg in args])

            # Get a unique name for the compiled function
            cname = "__{}_{}".format(
                func.__name__, hex(hash(arg_types) % ((sys.maxsize + 1) * 2))
            )

            # Compile the function if needed & cache it
            if not hasattr(instance, cname):

                dummy_args = [arg_type() for arg_type in arg_types]

                # Compile the function
                with CompileLogMessage(func.__name__):
                    with change_flags(compute_test_value="off"):
                        compiled_func = theano.function(
                            [*dummy_args],
                            func(instance, *dummy_args),
                            on_unused_input="ignore",
                            profile=config.profile,
                            mode=config.mode,
                        )
                    setattr(instance, cname, compiled_func)

            # d3v.d3viz(getattr(instance, cname), func.__name__+".html")

            # Return the compiled version
            return getattr(instance, cname)(*args)

    return wrapper


def clear_cache(instance, func):
    """
    Clear the compiled function cache for method `func` of a class
    instance `instance`.

    """
    basename = "__{}_".format(func.__name__)
    for key in list(instance.__dict__.keys()):
        if key.startswith(basename):
            delattr(instance, key)
