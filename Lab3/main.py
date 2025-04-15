import numpy as np
from scipy.optimize import linprog


def get_duality(A, b, c):
    # Коректне формування двоїстої задачі:
    # -b для мінімізації (оскільки linprog може лише мінімізувати)
    # A.T з оберненою нерівністю (≥ замість ≤)
    return -A.T, -c, -b


def print_problem(A, b, c):
    print("A:\n", A)
    print("b:\n", b)
    print("c:\n", c)


if __name__ == "__main__":
    c = np.array([1, 2])
    A = np.array(
        [
            [5, -2],
            [-1, 2],
            [1, 1],
        ]
    )
    b = np.array([4, 4, 4])

    # Вирішення початкової задачі
    bounds = [(0, None)] * len(c)
    result = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
    print("result:\n", result)

    # Отримання та вирішення двоїстої задачі
    A_d, b_d, c_d = get_duality(A, b, c)
    # Для двоїстої задачі ми використовуємо A_ub=-A.T та b_ub=-c,
    # оскільки linprog працює з обмеженнями виду A_ub*x <= b_ub
    d_bounds = [(0, None)] * len(c_d)
    d_result = linprog(-c_d, A_ub=-A_d, b_ub=-b_d, bounds=d_bounds, method="highs")

    # Перетворимо результат назад, щоб отримати правильне значення цільової функції
    d_result.fun = -d_result.fun

    print("duality_result:\n", d_result)
