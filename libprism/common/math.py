from numpy import asarray, sum, log, exp, rollaxis

def logsumexp(a, axis=None):
    """Compute the log of the sum of exponentials of input elements.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int, optional
        Axis over which the sum is taken. By default `axis` is None,
        and all elements are summed.

    Returns
    -------
    res : ndarray
        The result, ``np.log(np.sum(np.exp(a)))`` calculated in a numerically
        more stable way.

    See Also
    --------
    numpy.logaddexp, numpy.logaddexp2

    Notes
    -----
    Numpy has a logaddexp function which is very similar to `logsumexp`, but
    only handles two arguments. `logaddexp.reduce` is similar to this
    function, but may be less stable.

    Examples
    --------
    >>> from scipy.misc import logsumexp
    >>> a = np.arange(10)
    >>> np.log(np.sum(np.exp(a)))
    9.4586297444267107
    >>> logsumexp(a)
    9.4586297444267107

    """
    a = asarray(a)
    if axis is None:
        a = a.ravel()
    else:
        a = rollaxis(a, axis)
    a_max = a.max(axis=0)
    s = sum(exp(a - a_max), axis=0)
    # keyboard()
    if a.ndim == 1 and s == 0:
        return float("-inf")
    else:
        out = log(s)
    out += a_max
    return out#!/usr/bin/python
