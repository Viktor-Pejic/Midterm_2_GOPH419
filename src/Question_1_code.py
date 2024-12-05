





def cubic_spline(x):
    x = np.array(x, dtype=float)
    if np.any(x < xd[0]) or np.any(x > xd[-1]):
        raise ValueError(f"Input value(s) out of bounds: [{xd[0]}, {xd[-1]}].")
    n = len(xd) - 1
    div_dif_2 = np.diff(div_dif_1)

    rhs = np.zeros(n + 1)
    rhs[1:-1] = 3 * div_dif_2

    # Coefficient matrix
    A = np.zeros((n + 1, n + 1))
    A[0, 0], A[0, 1], A[0, 2] = delta_x[1], (delta_x[0] + delta_x[1]), -delta_x[0]
    A[-1, -3], A[-1, -2], A[-1, -1] = delta_x[-1], (delta_x[-2] + delta_x[-1]), -delta_x[-2]

    for i in range(1, n):
        A[i, i - 1] = delta_x[i - 1]
        A[i, i] = 2 * (delta_x[i - 1] + delta_x[i])
        A[i, i + 1] = delta_x[i]

    c = np.linalg.solve(A, rhs)

    d = np.diff(c) / (3 * delta_x)
    b = div_dif_1 - delta_x * (c[1:] + 2 * c[:-1]) / 3
    a = yd[:-1]

    # Evaluate spline
    i = np.searchsorted(xd, x) - 1
    i = np.clip(i, 0, n - 1)
    dx = x - xd[i]

    y = a[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
    return y