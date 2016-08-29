#! /usr/bin/env python

import sys
from ctypes import *


class Alignment(Structure):
    _fields_ = [
            ("seq1_aln" , c_char_p),
            ("seq2_aln" , c_char_p),
            ("aln_size",  c_int),
            ("score",     c_int)
        ]

BNW = CDLL('bnw.so')
BNW.align.argtypes = [ POINTER(c_char), c_int, POINTER(c_char), c_int, 
                   c_int, c_int, c_int, c_int, c_int ]

BNW.align.restype = POINTER(Alignment)
BNW.free_alignment.argtypes = [POINTER(Alignment)]



