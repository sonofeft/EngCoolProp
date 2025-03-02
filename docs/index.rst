
.. EngCoolProp documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../README.rst


EngCoolProp Interfaces The Coolprop Project with English Units
==============================================================

`CoolProp <http://www.coolprop.org/dev/index.html>`_ uses SI units for calculating
fluid properties. For those of us in industries that rely on English/Engineering
units, EngCoolProp provides an interface.

See the Code at: `<https://github.com/sonofeft/EngCoolProp>`_

See the Docs at: `<http://engcoolprop.readthedocs.org/en/latest/>`_

See PyPI page at:`<https://pypi.python.org/pypi/engcoolprop>`_

See the `CoolProp <http://www.coolprop.org/dev/index.html>`_ project at:
`<http://www.coolprop.org/dev/index.html>`_

EngCoolProp uses units of primarily inch, lbm, lbf, sec, BTU (some use of ft and hour).::

    The following are the default units for each property.

    T    = Temperature = degR
    Tc   = Critical Temperature = degR 
    Tnbp = Normal Boiling Point = degR
    P    = Pressure = psia
    Pc   = Critical Pressure = psia
    D    = Density = lbm/cu ft
    rho  = Density = lbm/cu inch
    Dc   = Critical Density = lbm/cu ft
    E    = Internal Energy = BTU/lbm
    H    = Enthalpy = BTU/lbm
    S    = Entropy = BTU/lbm degR
    Cv   = Heat Capacity (const. V) = BTU/lbm degR
    Cp   = Heat Capacity (const. P) = BTU/lbm degR
    g    = Ratio of Specific Heats (Cp/Cv) = (-)
    A    = Sonic Velocity = ft/sec
    V    = Viscosity = 1.0E5 * lb/ft-sec
    C    = Thermal Conductivity = BTU/ft-hr-R
    MW   = Molecular Weight = lbm/lbmmole
    Q    = Quality (fraction gas) = (-)
    Z    = Compressibility() = (-)



EngCoolProp
===========

Contents:

.. toctree::
   :maxdepth: 2

   quickstart
   examples
   incompressible
   incomp_solns
   functions
   incomp_fluid_funcs
   incomp_soln_funcs
   copyright
   authors
   history


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


