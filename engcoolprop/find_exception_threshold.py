import traceback


def find_exception_limit(f, tolerance=1e-6, lower_bound=0, upper_bound=10, 
                         no_exception='lower_bound', show_warnings=True):
    """
    Finds the value of x where f(x) does NOT throw an exception within the given tolerance.

    Parameters:
    f (function): The function to test.
    tolerance (float): The tolerance within which to find the value of x.
    lower_bound (float): The lower bound of the search interval.
    upper_bound (float): The upper bound of the search interval.

    no_exception(string): default="lower_bound": If no exceptions found return either lower_bound or upper_bound

    Returns:
    float: The value of x where f(x) does NOT throw an exception
    """
    tb_str = ''
    try:
        f( lower_bound )
        bad_lower = False
    except:
        bad_lower = True
        tb_str = traceback.format_exc()
        
    try:
        f( upper_bound )
        bad_upper = False
    except:
        bad_upper = True
        tb_str += '\n\n' + traceback.format_exc()
        
    # Exception boundary not included in bounds
    if bad_lower and bad_upper:
        if show_warnings:
            print( '*** Both upper and lower limits throw Exception in find_exception_limit***' )
            print( '    TRACEBACK:', tb_str)
        return None
        
    # All good
    if not bad_lower and not bad_upper:
        # print( '***No Exception found in find_exception_limit***' )

        # return either lower_bound or upper_bound if no exception found
        if no_exception=='lower_bound':
            return lower_bound
        else:
            return upper_bound
        
        
    while upper_bound - lower_bound > tolerance:
        mid = (lower_bound + upper_bound) / 2
        try:
            f(mid)
            if bad_lower:
                upper_bound = mid
            else:
                lower_bound = mid
        except Exception:
            if bad_lower:
                lower_bound = mid
            else:
                upper_bound = mid

    # Check and return the bound that does NOT throw an exception
    if bad_lower:
        return upper_bound
    else:
        return lower_bound
    #try:
    #    f(lower_bound)
    #    return lower_bound
    #except Exception:
    #    return upper_bound

if __name__ == "__main__":
    def f(x):
        # Dummy implementation of the function that may throw an exception.
        # Replace this with the actual function implementation as needed.
        if x >= 5:
            raise ValueError("x is too large")
        return x

    tolerance = 1e-6
    exception_value = find_exception_limit(f, tolerance, show_warnings=True)
    print(f"The value of x where f(x) throws an exception is approximately: {exception_value}")
