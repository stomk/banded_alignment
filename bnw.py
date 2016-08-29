import os
import sys
from ctypes import *

class Alignment(Structure):
    _fields_ = [
        ("seq1_aln" , c_char_p),
        ("seq2_aln" , c_char_p),
        ("aln_size",  c_int),
        ("score",     c_int)
        ]

dll_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib', 'bnw.so')
bnw = CDLL(dll_path)

bnw.align.argtypes = [ POINTER(c_char), c_int, POINTER(c_char), c_int, 
                   c_int, c_int, c_int, c_int, c_int ]

bnw.align.restype = POINTER(Alignment)
bnw.free_alignment.argtypes = [POINTER(Alignment)]

