import numpy as np

def coeffMatrix():
    CO2 = np.loadtxt("C:/Users/Viktor/repos/Midterm_2_GOPH419/data/CO2_data.txt", 'float')
    xd = np.array(CO2[51:62,0], dtype='float') #Year on the x-axis
    yd = np.array(CO2[51:62,1], dtype='float') #CO2 concentration on the y-axis

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

    return A

def main():
    A = coeffMatrix()
    print('\nCoefficient Matrix:\n')
    print(A)

if __name__ == '__main__':
    main()