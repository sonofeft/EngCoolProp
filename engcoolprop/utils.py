
import re
import CoolProp.CoolProp as CP

incomp_pure_solnL = CP.get_global_param_string('incompressible_list_solution').split(',')

incomp_pure_fluidL = CP.get_global_param_string('incompressible_list_pure').split(',')


def parse_coolprop_mixture(mixture_name):
    """
    Parses a CoolProp mixture name to extract the base name and mass fraction percentage.

    Args:
        mixture_name (str): The CoolProp mixture name (e.g., "LiBr[0.23]" or "LiBr-23%").

    Returns:
        tuple: A tuple containing the base name (str) and mass fraction percentage (float).
               Returns (None, None) if parsing fails.
    """
    if mixture_name is None or not isinstance(mixture_name, str):
        return None, None
    
    # Pattern for just the base name, mass fraction is 100%
    if mixture_name in incomp_pure_solnL:
        return mixture_name, 100.0 # for name w/o fraction_mass, assume 100%

    # Pattern for [fraction] notation
    match_bracket = re.match(r"([a-zA-Z0-9]+)\[([0-9.]+)\]", mixture_name)
    if match_bracket:
        base_name = match_bracket.group(1)
        mass_fraction_percentage = float(match_bracket.group(2)) * 100.0
        return base_name, mass_fraction_percentage

    # Pattern for -fraction% notation
    match_percent = re.match(r"([a-zA-Z0-9]+)-([0-9.]+)\%", mixture_name)
    if match_percent:
        base_name = match_percent.group(1)
        mass_fraction_percentage = float(match_percent.group(2))
        return base_name, mass_fraction_percentage

    return None, None  # Parsing failed

    # # Example usage:
    # print(parse_coolprop_mixture("LiBr[0.23]"))
    # print(parse_coolprop_mixture("LiBr-23%"))
    # print(parse_coolprop_mixture("Water"))
    # print(parse_coolprop_mixture("NaCl[0.1]"))
    # print(parse_coolprop_mixture("SomeInvalidString"))

def print_avoid_valerr():
    print( '.'*22, 'To avoid ValueError set "auto_fix_value_errors" to True', '.'*22 )

def format_float(float_val, output_len=9, sig_digits=4):
    """
    Formats a floating point number to a string of a specified length and significant digits.

    Parameters:
    float_val (float): The floating point number to format.
    output_len (int): The required length of the output string.
    sig_digits (int): The number of significant digits of float_val in the output string.

    Returns:
    str: The formatted string.
    """

    #  using %g format and cleaning up long exponents if "e" format
    format_str = f"%.{sig_digits}g"
    result = format_str % float_val
    result = result.replace('e+0', 'e').replace('e-0', 'e-').replace('e+', 'e')
    
    return result.ljust(output_len)

class Same_g_len:
    def __init__(self, obj, properties):
        """
        Initialize Same_g_len with an object and a list of property names.

        Parameters:
        obj (object): The object from which to get the property values.
        properties (list): A list of strings that are the names of the object's properties.
        
        Returns: properties padded with spaces to _max_length
        """
        self._max_length = 0
        self._properties = {}

        for prop in properties:
            if hasattr(obj, prop):
                value = getattr(obj, prop)
                try:
                    formatted_value = f"{value:g}"
                except:
                    formatted_value = str( value )
                self._properties[prop] = formatted_value
                self._max_length = max(self._max_length, len(formatted_value))

    def __getattr__(self, name):
        """
        Override __getattr__ to return right justified string of the maximum length.

        Parameters:
        name (str): The name of the attribute to access.

        Returns:
        str: The right justified string of the maximum length.
        """
        if name in self._properties:
            return self._properties[name].rjust(self._max_length)
        else:
            return '?'*self._max_length

def in_between(value, val1, val2):
    """
    Check if a value is between val1 and val2 inclusive, regardless of the order of val1 and val2.

    Parameters:
    value (float): The value to check.
    val1 (float): One end of the range.
    val2 (float): The other end of the range.

    Returns:
    bool: True if value is between val1 and val2 inclusive, False otherwise.
    """
    min_val = min(val1, val2)
    max_val = max(val1, val2)
    return min_val <= value <= max_val

def opposite_signs(x, y):
    """
    Check if two numbers have opposite signs.
    Zero is considered to have the opposite sign to a negative number.

    Parameters:
    x (float): The first number.
    y (float): The second number.

    Returns:
    bool: True if x and y have opposite signs, False otherwise.
    """
    return (x == 0 and y != 0) or (y == 0 and x != 0) or (x < 0 < y) or (y < 0 < x)

if __name__ == "__main__":

    # Example usage
    print( '... format_float ...')
    print('"%s"'%format_float(12345678.9))  # Output: "1.235e7  "
    print('"%s"'%format_float(0.00000123))  # Output: "1.23e-6  "
    print('"%s"'%format_float(-9876543.21)) # Output: "-9.877e6 "
    print('"%s"'%format_float(123.456))     # Output: "123.5    "
    print('"%s"'%format_float(0.123456))    # Output: "0.1235   "


    # Example usage
    print()
    print( '... in_between ...')
    print(in_between(5, 1, 10))  # Output: True
    print(in_between(5, 10, 1))  # Output: True
    print(in_between(1, 10, 1))  # Output: True
    print(in_between(10, 1, 10)) # Output: True
    print(in_between(0, 1, 10))  # Output: False
    print(in_between(15, 1, 10)) # Output: False

    # Example usage
    print()
    print( '... opposite_signs ...')
    # Example usage with expected answers
    def test_opposite_signs():
        assert opposite_signs(5, -10) == True, "Test case 1 failed"
        assert opposite_signs(-5, 10) == True, "Test case 2 failed"
        assert opposite_signs(0, -10) == True, "Test case 3 failed"
        assert opposite_signs(-10, 0) == True, "Test case 4 failed"
        assert opposite_signs(10, 0) == True, "Test case 5 failed"
        assert opposite_signs(5, 10)   == False, "Test case 6 failed"
        assert opposite_signs(-5, -10) == False, "Test case 7 failed"
        assert opposite_signs(0, 0)    == False, "Test case 8 failed"
        print("All test cases passed")

    # Run the tests
    test_opposite_signs()

