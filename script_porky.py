from pymol.cgo import *
from pymol import cmd

cmd.load("avg_1hvr_apo.pdb")
cmd.load("disp_1hvr.pdb.pdb")
cmd.load("modevectors.py")
modevectors("avg_1hvr_apo", "disp_1hvr.pdb", outname="modevectors", head=1.0, tail=0.3, headrgb = "1.0, 1.0, 0.0", tailrgb = "1.0, 1.0, 0.0") 