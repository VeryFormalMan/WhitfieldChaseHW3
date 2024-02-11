# Chase Whitfield
# MAE 3403
# Homework 3 Part B
# 02/12/2024

import math

def t_distribution_probability(degrees_of_freedom, z_value):
    """
    Calculate probability using t-distribution formula.
    T-distribution is probability distribution used in hypothesis
    testing when sample size is small and population standard deviation is unknown.
    :param degrees_of_freedom: Degrees of freedom of t-distribution.
    Represents number of independent observations used to estimate a parameter.
    :param z_value: The value of z for which we want to calculate the probability.
    :return: float: Calculated probability based on the t-distribution.
    """
    numerator = math.gamma((degrees_of_freedom + 1) / 2)
    denominator = math.sqrt(degrees_of_freedom * math.pi) * math.gamma(degrees_of_freedom / 2)
    coefficient = numerator / denominator
    exponent = -((degrees_of_freedom + 1) / 2)
    probability = coefficient * ((1 + (z_value ** 2) / degrees_of_freedom) ** exponent)
    return probability

# Used GPT to help create docfile and to help with lines 18-20
def main():
    """
    Main function.
    :return: Displays probability.
    """
    degrees_of_freedom = int(input("Enter the degrees of freedom: "))
    z_value = float(input("Enter the value of z: "))
    probability = t_distribution_probability(degrees_of_freedom, z_value)
    print("Probability:", probability)

if __name__ == "__main__":
    main()
