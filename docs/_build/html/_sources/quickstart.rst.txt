
.. quickstart

QuickStart
==========

Install EngCoolProp
-------------------

The easiest way to install EngCoolProp is::

    pip install engcoolprop
    
        OR on Linux
    sudo pip install engcoolprop
        OR perhaps
    pip install --user engcoolprop

In case of error, see :ref:`internal_pip_error`

.. _internal_source_install:

Installation From Source
------------------------

Much less common, but if installing from source, then
the best way to install engcoolprop is still ``pip``.

After navigating to the directory holding EngCoolProp source code, do the following::

    cd full/path/to/engcoolprop
    pip install -e .
    
        OR on Linux
    sudo pip install -e .
        OR perhaps
    pip install --user -e .
    
This will execute the local ``setup.py`` file and insure that the pip-specific commands in ``setup.py`` are run.

Running EngCoolProp
-------------------

A simple usage of EngCoolProp is shown here::

    # import the package
    from engcoolprop.ec_fluid import EC_Fluid

    # create a EC_Fluid object (possibly intialize with T and P)
    ec = EC_Fluid(symbol="N2", T=530.0,P=100.0 ) # T=degR, P=psia
    
    # set to desired properties
    ec.setProps(T=500., D=0.1)

    # print a summary of the properties
    ec.printProps()


.. _internal_pip_error:

pip Error Messages
------------------

If you get an error message that ``pip`` is not found, see `<https://pip.pypa.io/en/latest/installation/>`_ for full description of ``pip`` installation.

There might be issues with ``pip`` failing on Linux with a message like::


    InsecurePlatformWarning
            or    
    Cannot fetch index base URL https://pypi.python.org/simple/

Certain Python platforms (specifically, versions of Python earlier than 2.7.9) have the InsecurePlatformWarning. 
If you encounter this warning, it is strongly recommended you upgrade to a newer Python version, or that you use pyOpenSSL.    

Also ``pip`` may be mis-configured and point to the wrong PyPI repository.
You need to fix this global problem with ``pip`` just to make python usable on your system.


If you give up on upgrading python or fixing ``pip``, 
you might also try downloading the engcoolprop source package 
(and all dependency source packages)
from PyPI and installing from source as shown above at :ref:`internal_source_install`


