
import numpy as np
import matplotlib.pyplot as plt



def cubic_spline(x, xd, yd):
    x = np.array(x, dtype=float)
    if np.any(x < xd[0]) or np.any(x > xd[-1]):
        raise ValueError(f"Input value(s) out of bounds: [{xd[0]}, {xd[-1]}].")
    n = len(xd) - 1

    delta_x = np.diff(xd)
    delta_y = np.diff(yd)
    div_dif_1 = delta_y / delta_x
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

def main():
    CO2 = np.loadtxt("C:/Users/Viktor/repos/Midterm_2_GOPH419/data/CO2_data.txt", 'float')
    xd = np.array(CO2[51:62,0], dtype='float') #Year on the x-axis
    yd = np.array(CO2[51:62,1], dtype='float') #CO2 concentration on the y-axis

    x = np.linspace(2010,2020,10)

    spline = cubic_spline(x, xd, yd)

    plt.plot(x, spline)
    plt.scatter(xd, yd, color='black')
    plt.xlabel('Years')
    plt.ylabel('CO2 Concentration (mean)')
    plt.title('CO2 Concentration Trend from 2010-2020', fontweight='bold')
    plt.show()

if __name__ == '__main__':
    main()
