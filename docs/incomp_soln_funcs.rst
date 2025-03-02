
.. incomp_soln_funcs


.. _link_ec_incomp_soln_functions:

EC_Incomp_Soln Functions
========================

Create a EC_Incomp_Soln object like this::
    
    from engcoolprop.ec_incomp_soln import EC_Incomp_Soln

    # Create incompressible object. (without specifying state point)
    ec_soln = EC_Incomp_Soln(symbol="MEG-30%")

If the mass fraction is not known, simply use the base name to get an average
supported mass fraction like this::

    ec_soln = EC_Incomp_Soln(symbol="MEG")

**EC_Incomp_Soln Description**

.. automodule:: engcoolprop.ec_incomp_soln

**EC_Incomp_Soln Methods**

.. autoclass:: EC_Incomp_Soln
   :members:
