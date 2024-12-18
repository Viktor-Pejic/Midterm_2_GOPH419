import numpy as np
import matplotlib.pyplot as plt

CO2 = np.loadtxt("C:/Users/Viktor/repos/Midterm_2_GOPH419/data/CO2_data.txt", 'float')
xd = np.array(CO2[:, 0], dtype='float')  # Year 2010-2020 on the x-axis
yd = np.array(CO2[:, 1], dtype='float')  # CO2 concentration on the y-axis

def second_order(xd, yd):
    h = xd[1] - xd[0]
    f_2prime_2 = np.zeros(len(yd))
    f_2prime_2[0], f_2prime_2[-1] = np.nan, np.nan
    for i in range(1, len(yd) - 1):
        f_2prime_2[i] = (yd[i + 1] - 2*yd[i] + yd[i - 1]) / (h**2)
    return f_2prime_2

def fourth_order(xd, yd):
    h = xd[1] - xd[0]
    f_2prime_4 = np.zeros(len(yd))
    f_2prime_4[0], f_2prime_4[1], f_2prime_4[-1], f_2prime_4[-2] = np.nan, np.nan, np.nan, np.nan
    for i in range(2, len(yd) - 2):
        f_2prime_4[i] = (-yd[i + 2] + 16*yd[i + 1] - 30*yd[i] + 16*yd[i - 1] - yd[i - 2]) / (12*(h**2))
    return f_2prime_4

def main(xd, yd):
    order_2 = second_order(xd, yd)
    order_4 = fourth_order(xd, yd)


    plt.plot(xd, order_2, label='second_order')
    plt.plot(xd, order_4, label='fourth_order')
    plt.xlabel('Year')
    plt.ylabel('CO2 Concentration / Rate of Change (ppm/yr^2)')
    plt.title('CO2 Concentration and its First Derivative')
    plt.legend()
    plt.grid()
    plt.savefig("C:/Users/Viktor/repos/Midterm_2_GOPH419/figures/Second derivative 2nd and 4th order 1959-2022.png")

if __name__ == '__main__':
    main(xd, yd)



