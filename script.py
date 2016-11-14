from pymol.cgo import *
from pymol import cmd

cmd.load("avg_1hvr_apo.pdb")
cmd.load("out.pdb")
cmd.load("modevectors.py")
modevectors("avg_1hvr_apo", "out")
