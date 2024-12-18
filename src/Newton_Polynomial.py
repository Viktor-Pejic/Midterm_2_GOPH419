import numpy as np


def newton_dd(xd, yd):
    n = len(xd)
    z = np.zeros((n, n))
    z[:, 0] = yd

    coeff = np.zeros(n)
    coeff[0] = yd[0]

    # Fill the divided difference table
    for j in range(1, n):
        for i in range(n - j):
            z[i, j] = (z[i + 1, j - 1] - z[i, j - 1]) / (xd[i + j] - xd[i])
        coeff[j] = z[0, j]

    return coeff


def newton_polynomial(x, xd, coeff):
    n = len(xd)
    result = coeff[0]

    for i in range(1, n):
        term = coeff[i]
        for j in range(i):
            term *= (x - xd[j])
        result += term

    return result


def estimate_concentration(xd, yd, target_x, tolerance=1e-4):
    order = 2
    previous_value = None
    while True:
        coeff = newton_dd(xd[:order], yd[:order])  # Compute coefficients for the current order
        estimated_value = newton_polynomial(target_x, xd[:order], coeff)  # Estimate value
        if previous_value is not None and abs(estimated_value - previous_value) < tolerance:
            break
        previous_value = estimated_value
        order += 1

    return estimated_value, order

def find_sufficient_degree(xd, yd, target_x, tolerance=0.0005):
    previous_value = None
    for degree in range(2, len(xd)):  # Start from degree 2
        coeff = newton_dd(xd[:degree + 1], yd[:degree + 1])  # Compute coefficients
        current_value = newton_polynomial(target_x, xd[:degree + 1], coeff)  # Estimate

        if previous_value is not None and abs(current_value - previous_value) < tolerance:
            return current_value, degree  # Return the value and degree
        previous_value = current_value

    return current_value, degree  # Return the highest degree if tolerance is not met


def main():

    CO2 = np.loadtxt("C:/Users/Viktor/repos/Midterm_2_GOPH419/data/CO2_data.txt", 'float')
    xd = np.array(CO2[51:62, 0], dtype='float')  # Year 2010-2020 on the x-axis
    yd = np.array(CO2[51:62, 1], dtype='float')  # CO2 concentration on the y-axis


    target_x = 2015.25


    estimated_value, order = estimate_concentration(xd, yd, target_x)

    print(f"Estimated CO2 concentration for March 2015: {estimated_value:.4f}")
    print(f"Polynomial order required: {order}\n")

    result, degree = find_sufficient_degree(xd, yd, target_x)
    print(f"Estimated CO2 concentration up to 6 significant digits: {result:.3f}")
    print(f"Required polynomial degree: {degree}")

if __name__ == '__main__':
    main()