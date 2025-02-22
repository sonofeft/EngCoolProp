import numpy as np


class RootResults:
    def __init__(self, root, iterations, function_calls, flag, converged):
        self.root = root
        self.iterations = iterations
        self.function_calls = function_calls
        self.flag = flag
        self.converged = converged

    def __repr__(self):
        return (f"RootResults(root={self.root}, iterations={self.iterations}, "
                f"function_calls={self.function_calls}, flag={self.flag}, "
                f"converged={self.converged})")


# brentq method take from scipy
# (modified slightly for initial lower and upper bounds within xtol)
def brentq(f, a, b, args=(), xtol=1e-12, rtol=4*np.finfo(float).eps, maxiter=100, **kwargs):
    """
    Brent's method for root finding.
    
    Parameters
    ----------
    f : function
        Python function returning a scalar. The root to be found is such that f(x) = 0.
    a : scalar
        Lower bound of the bracketing interval [a,b].
    b : scalar
        Upper bound of the bracketing interval [a,b].
    xtol : scalar, optional
        Absolute error tolerance. Default is 1e-12.
    rtol : scalar, optional
        Relative error tolerance. Default is 4 * machine epsilon for float.
    maxiter : int, optional
        Maximum number of iterations. Default is 100.
    args : tuple, optional
        Extra arguments for the function `f`.
    
    """
    iters = 0
    funcalls = 1
    fa = f(a, *args)
    if abs(fa) < xtol:
        return RootResults(root=a, iterations=iters, function_calls=funcalls, flag=0, converged=True)
    
    funcalls = 2
    fb = f(b, *args)
    if abs(fb) < xtol:
        return RootResults(root=b, iterations=iters, function_calls=funcalls, flag=0, converged=True)

    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have different signs")

    c, fc = a, fa
    d = e = b - a

    for iters in range(maxiter):
        if fb * fc > 0:
            c, fc, d, e = a, fa, b - a, b - a

        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb

        tol = 2 * xtol * max(abs(b), 1.0) + 0.5 * rtol
        m = 0.5 * (c - b)

        if abs(m) <= tol or fb == 0:
            return RootResults(root=b, iterations=iters, function_calls=funcalls, flag=0, converged=True)

        if abs(e) >= tol and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2 * m * s
                q = 1 - s
            else:
                q = fa / fc
                r = fb / fc
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)

            if p > 0:
                q = -q
            p = abs(p)

            if 2 * p < min(3 * m * q - abs(tol * q), abs(e * q)):
                e, d = d, p / q
            else:
                d, e = m, m
        else:
            d, e = m, m

        a, fa = b, fb
        if abs(d) > tol:
            b += d
        else:
            b += np.sign(m) * tol

        fb = f(b, *args)
        funcalls += 1

    return RootResults(root=b, iterations=iters, function_calls=funcalls, flag=1, converged=False)
        

# Example usage
if __name__ == "__main__":

    def root_function(x):
        """
        Define the function for which we want to find the root.
        """
        # Example function: f(x) = x^2 - 2
        return x**2 - 2

    def find_root_brentq(func, a, b, tol=1e-5):
        """
        Find the root of the given function using the Brentq method.
        
        Parameters:
        func (callable): The function for which to find the root.
        a (float): The lower bound of the interval.
        b (float): The upper bound of the interval.
        tol (float): The tolerance for the root-finding algorithm. Default is 1e-5.
        
        Returns:
        float: The root of the function within the given interval.
        """
        sol = brentq(root_function, 0, 2)
        print(sol)
        return sol.root


    a = 0  # Lower bound of the interval
    b = 2  # Upper bound of the interval
    root = find_root_brentq(root_function, a, b)
    print(f"The root of the function in the interval [{a}, {b}] is: {root}")        
