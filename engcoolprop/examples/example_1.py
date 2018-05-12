"""
This example demonstrates the proper use of project: engcoolprop
"""
from __future__ import print_function
import sys
import os

sys.path.insert(0, os.path.abspath("../../"))  # needed to find engcoolprop development version

from engcoolprop.ec_fluid import EC_Fluid

ec = EC_Fluid(symbol="N2", T=530.0,P=100.0 ) # T=degR, P=psia
ec.setProps(T=500., D=0.1)
ec.printTPD()
print('')

ec.printCriticalProps()
print('')

ec.printProps()

