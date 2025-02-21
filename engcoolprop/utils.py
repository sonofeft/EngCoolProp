

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

