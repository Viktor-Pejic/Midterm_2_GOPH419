import numpy as np
from src.Coefficient_Matrix import coeffMatrix

def GaussSeidel(A, xd, yd, x0=None, tol=1e-8, max_iter=100):
    n = len(xd) - 1

    delta_x = np.diff(xd)
    delta_y = np.diff(yd)
    div_dif_1 = delta_y / delta_x
    div_dif_2 = np.diff(div_dif_1)

    rhs = np.zeros(n + 1)
    rhs[1:-1] = 3 * div_dif_2

    A = np.array(A, dtype=float)
    b = np.array(rhs, dtype=float)

    # Initialize x0 if None
    if x0 is None:
        x = np.zeros_like(b, dtype=float)
    else:
        x0 = np.array(x0, dtype=float)
        if x0.shape != b.shape:
            raise ValueError("Initial guess x0 must have the same shape as b.")
        x = x0.copy()

    n = A.shape[0]

    for iteration in range(max_iter):
        x_old = x.copy()
        for i in range(n):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i+1:], x[i+1:])
            x[i] = (b[i] - sum1 - sum2) / A[i, i]

        # Check for convergence
        if np.max(np.abs(x - x_old)) / np.max(np.abs(x)) < tol:
            expected_sol = np.linalg.solve(A, rhs)
            for i in range(len(x)):
                diff = abs(x[i] - expected_sol[i])
                if diff <= 1e-03:
                    print("\nGauss-Seidel solution satisfies the original equation\n")
                    return x
                else:
                    print("\nGauss-Seidel solution does not satisfy the original equation\n")


    # If max_iter reached without convergence
    raise RuntimeWarning(f"Solution did not converge within {max_iter} iterations.")

def main():
    CO2 = np.loadtxt("C:/Users/Viktor/repos/Midterm_2_GOPH419/data/CO2_data.txt", 'float')
    xd = np.array(CO2[51:62, 0], dtype='float')  # Year 2010-2020 on the x-axis
    yd = np.array(CO2[51:62, 1], dtype='float')  # CO2 concentration on the y-axis

    A = coeffMatrix()

    solution = GaussSeidel(A, xd, yd, x0=None, tol=1e-8, max_iter=100)

    print(solution)


if __name__ == "__main__":
    main()
