
import re
import CoolProp.CoolProp as CP

incomp_pure_fluidL = CP.get_global_param_string('incompressible_list_solution').split(',')

def parse_coolprop_mixture(mixture_name):
    """
    Parses the name of a CoolProp mixture and returns the base name and mass percentage.

    Args:
    - mixture_name (str): The name of the CoolProp mixture (e.g., "LiBr[0.23]", "LiBr-23%", or "MAM2-23.6%").

    Returns:
    - tuple: A tuple containing the base name (str) and the mass percentage (float).
    """
    if mixture_name in incomp_pure_fluidL:
        return mixture_name, 100 # for name w/o fraction_mass, assume 100%

    # Regular expression patterns to match the different formats
    pattern_brackets = re.compile(r"(?P<base_name>[A-Za-z0-9]+)\[(?P<mass_fraction>\d*\.\d+|\d+)\]")
    pattern_percentage = re.compile(r"(?P<base_name>[A-Za-z0-9]+)-(?P<mass_fraction>\d*\.\d+|\d+)%")

    match_brackets = pattern_brackets.match(mixture_name)
    match_percentage = pattern_percentage.match(mixture_name)

    if match_brackets:
        base_name = match_brackets.group("base_name")
        mass_fraction = float(match_brackets.group("mass_fraction"))
        percent_mass = mass_fraction * 100
    elif match_percentage:
        base_name = match_percentage.group("base_name")
        percent_mass = float(match_percentage.group("mass_fraction"))
    else:
        raise ValueError("Invalid mixture name format")

    return base_name, percent_mass


    # # Example usage
    # mixture_name_1 = "LiBr[0.23]"
    # mixture_name_2 = "LiBr-23%"

    # base_name, percent_mass = parse_coolprop_mixture(mixture_name_1)
    # base_name, percent_mass = parse_coolprop_mixture(mixture_name_2)

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

