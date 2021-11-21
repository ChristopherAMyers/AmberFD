def AssertEqual(val1, val2, tol):
    denom = max(abs(val1), abs(val2))
    diff = abs(val1 - val2)/denom
    assert(diff < tol)