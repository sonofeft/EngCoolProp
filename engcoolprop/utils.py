

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


if __name__ == "__main__":

    # Example usage
    print('"%s"'%format_float(12345678.9))  # Output: "1.235e7  "
    print('"%s"'%format_float(0.00000123))  # Output: "1.23e-6  "
    print('"%s"'%format_float(-9876543.21)) # Output: "-9.877e6 "
    print('"%s"'%format_float(123.456))     # Output: "123.5    "
    print('"%s"'%format_float(0.123456))    # Output: "0.1235   "




