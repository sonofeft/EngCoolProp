
import unittest
# import unittest2 as unittest # for versions of python < 2.7

"""
        Method                            Checks that
self.assertEqual(a, b)                      a == b   
self.assertNotEqual(a, b)                   a != b   
self.assertTrue(x)                          bool(x) is True  
self.assertFalse(x)                         bool(x) is False     
self.assertIs(a, b)                         a is b
self.assertIsNot(a, b)                      a is not b
self.assertIsNone(x)                        x is None 
self.assertIsNotNone(x)                     x is not None 
self.assertIn(a, b)                         a in b
self.assertNotIn(a, b)                      a not in b
self.assertIsInstance(a, b)                 isinstance(a, b)  
self.assertNotIsInstance(a, b)              not isinstance(a, b)  
self.assertAlmostEqual(a, b, places=5)      a within 5 decimal places of b
self.assertNotAlmostEqual(a, b, delta=0.1)  a is not within 0.1 of b
self.assertGreater(a, b)                    a is > b
self.assertGreaterEqual(a, b)               a is >= b
self.assertLess(a, b)                       a is < b
self.assertLessEqual(a, b)                  a is <= b

for expected exceptions, use:

with self.assertRaises(Exception):
    blah...blah...blah

with self.assertRaises(KeyError):
    blah...blah...blah

Test if __name__ == "__main__":


    def test__main__(self):
        # Change bottom of source file to call "dev_tests"
        
         def dev_tests():
            pass

         if __name__ == "__main__":
            dev_tests()
            
        # then test by calling <name>.dev_tests()


See:
      https://docs.python.org/2/library/unittest.html
         or
      https://docs.python.org/dev/library/unittest.html
for more assert options
"""

import sys, os

from engcoolprop.ec_humid_air import EC_Humid_Air, dev_tests
from engcoolprop.humid_air_params import *

here = os.path.abspath(os.path.dirname(__file__)) # Needed for py.test
up_one = os.path.split( here )[0]  # Needed to find models development version
if here not in sys.path[:2]:
    sys.path.insert(0, here)
if up_one not in sys.path[:2]:
    sys.path.insert(0, up_one)

class MyTest(unittest.TestCase):


    def test_should_always_pass_cleanly(self):
        """Should always pass cleanly."""
        pass

    def test_myclass_existence(self):
        """Check that myclass exists"""
        ec_ha = EC_Humid_Air()
        
        # See if the self.EC_Incomp_Fluid object exists
        self.assertIsInstance(ec_ha, EC_Humid_Air, msg=None)

    def test_print_IO(self):
        """test print IO"""

        ec_ha = EC_Humid_Air()
        ec_ha.print_input_params()
        ec_ha.print_output_params()
        ec_ha.printProps()
        ec_ha.printProps( eng_units=False )


    def test_bad_input_pair(self):
        """test bad_input_pair"""
        
        ec_ha = EC_Humid_Air( Tdb=530.0, RelHum=0.5 )

        ec_ha2 = EC_Humid_Air( HumRat=ec_ha.HumRat, RelHum=ec_ha.RelHum )
        self.assertEqual(  ec_ha2.Cond, float('inf') )
        

    def test_degF(self):
        """test degF"""
        ec_ha = EC_Humid_Air( Tdb=530.0, RelHum=0.5 )

        ec_ha2 = EC_Humid_Air(  TdbdegF=70.33, RelHum=0.5  )
        
        self.assertAlmostEqual(ec_ha2.Cond, ec_ha.Cond, places=5)

    def test_wrong_number_of_inputs(self):
        """test wrong number of inputs"""

        with self.assertRaises(ValueError):
            ec_ha = EC_Humid_Air( Tdb=530.0 )


    def test_omitted_pressure(self):
        """test omitted pressure"""
        
        with self.assertRaises(ValueError):
            ec_ha = EC_Humid_Air( Tdb=530.0, WetBulb=520, DewPoint=510 )

    def test_two_water_contents(self):
        """test two water contents"""
        
        ec_ha = EC_Humid_Air( Tdb=530.0, RelHum=0.5 )

        ec_ha2 = EC_Humid_Air(  RelHum=ec_ha.RelHum, DewPoint=ec_ha.DewPoint  )
        
        with self.assertRaises(AssertionError):
            self.assertAlmostEqual(ec_ha2.Cond, ec_ha.Cond, places=5)

    def test__main__(self):
        old_sys_argv = list(sys.argv)
        sys.argv = list(sys.argv)
        sys.argv.append('suppress_show')
        
        try:
            dev_tests()
        except:
            raise Exception('ERROR... failed in __main__ routine')
        finally:
            sys.argv = old_sys_argv
        

if __name__ == '__main__':
    # Can test just this file from command prompt
    #  or it can be part of test discovery from nose, unittest, pytest, etc.
    unittest.main()

