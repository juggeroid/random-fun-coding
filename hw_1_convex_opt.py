import numpy
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.axislines import SubplotZero
import cvxpy

# Variable constants (points, coordinates).
X_MIN, X_MAX = -1, 1
Y_MIN, Y_MAX = -1, 1

def generate_points(points, low_variance = True):
    x = np.random.uniform(X_MIN, X_MAX, size = points)
    if low_variance:
        # Scale the points so that they'll be distributed along the artificial polynomial.
        # y = 3x + 1 + epsilon
        errors = np.random.uniform(Y_MIN, Y_MAX, size = points)
        y = 1 + 3 * x + errors * np.random.normal(size = points)
    else:
        y = np.random.uniform(Y_MIN, Y_MAX, size = points)
    # Sort the points' coordinates by x in an increasing order.
    xs, ys = map(numpy.array, zip(*sorted(zip(x, y), reverse = False)))
    return xs, ys

def plot_points_helper(x, y, x_poly, y_poly, scale_plot = True):

    figure = plt.figure(figsize = (12, 8))
    ax = SubplotZero(figure, 111)
    figure.add_subplot(ax)

    for direction in ["xzero", "yzero"]:
        ax.axis[direction].set_axisline_style("-|>")
        ax.axis[direction].set_visible(True)

    for direction in ["left", "right", "bottom", "top"]:
        ax.axis[direction].set_visible(False)

    scale = np.mean(x[x > 0]) + np.mean(y[y > 0]) if scale_plot else 0
    plt.xlim([min(x) - scale, max(x) + scale])
    plt.ylim([min(y) - scale, max(y) + scale])

    plt.scatter(x, y)
    plt.plot(x_poly, y_poly, color = 'r')
    plt.show()

def evaluate_polynomial(coefficients):
    # Generate the span for which the polynomial will be constructed.
    x = np.linspace(X_MIN, X_MAX, num = 500)
    # Form the polynomial by element-wise multiplication, raise x to the power of p and multiply by the coefficient.
    y = coefficients @ [np.power(x, p) for p in range(len(coefficients))]
    return x, y

def polynomial_approximation(x, y, degree):
    # Degree is shifted so the addition has to be performed.
    degree += 1
    A = numpy.matrix([np.power(x, p) for p in range(degree)])
    # Transpose the matrix to acquire the correct dimension.
    A = A.T
    coefficients = cvxpy.Variable(degree)
    # Formulate the problem of approximation as described here.
    # https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf
    problem_formulation = cvxpy.sum_squares(A @ coefficients - y)
    problem_instance = cvxpy.Problem(cvxpy.Minimize(problem_formulation))
    error = problem_instance.solve(solver = cvxpy.SCS)
    values = coefficients.value
    return error, values

def get_all_approximation_errors(x, y):
    errors = []
    solutions = []
    for degree in range(len(y)):
        error, solution = polynomial_approximation(x, y, degree = degree)
        errors.append(error)
        solutions.append(solution)
    return np.array(errors), solutions

def dirty_test_function(points = 10, low_variance = False):
    x, y = generate_points(points, low_variance = low_variance)
    errors, solutions = get_all_approximation_errors(x, y)
    plt.plot(errors[1:])
    print(errors[1:])
    plot_points_helper(x, y, *evaluate_polynomial(solutions[0]))
    plot_points_helper(x, y, *evaluate_polynomial(solutions[1]))
    plot_points_helper(x, y, *evaluate_polynomial(solutions[2]))
    plot_points_helper(x, y, *evaluate_polynomial(solutions[-1]))

dirty_test_function(10, low_variance = False)
dirty_test_function(10, low_variance = True)