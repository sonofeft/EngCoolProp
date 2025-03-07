
.. humid_air_funcs


.. _link_ec_humid_air_functions:

EC_Humid_Air Functions
======================

Create a EC_Humid_Air object like this::
    
    from engcoolprop.ec_humid_air import EC_Humid_Air

    # Create incompressible object. (without specifying state point)
    ec_ha = EC_Humid_Air()
       ... OR ...
    ec_ha = EC_Humid_Air(T=536.4, P=14.6959, RelHum=0.5)

**EC_Humid_Air Description**

.. automodule:: engcoolprop.ec_humid_air

**EC_Humid_Air Methods**

.. autoclass:: EC_Humid_Air
   :members:
