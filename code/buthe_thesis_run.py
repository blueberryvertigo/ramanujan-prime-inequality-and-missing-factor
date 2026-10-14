"""Refereed-input Theorem A' and Theorem D (buthe_ramanujan.py) with Buthe's Theorem 1 as corrected in
Bhattacharjee's thesis (E2 of (2.15), E3 of (2.16)) instead of the log(3c)-only correction."""
import runpy, buthe_bound as BB
BB.E2MODE = 'thesis'
runpy.run_module('buthe_ramanujan', run_name='__main__')
