def AssertEqual(val1, val2, tol, print_result=False):
    denom = max(abs(val1), abs(val2))
    diff = abs(val1 - val2)/denom
    if print_result or (diff > tol):
        result_str = "val1={:.16e}; val2={:.16e}; diff={:.16e}; tol={:.16e}".format(val1, val2, diff, tol)
        if (diff > tol):
            raise Exception("Assertion Failed; " + result_str)
        else:
            print("Assertion: " + result_str)        
    #assert(diff < tol)