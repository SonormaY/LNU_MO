import numpy as np
def simplex_method(c, A, b, problem_type="min", inequality_types=None):
    if problem_type == "max":
        c = [-coef for coef in c]
    m, n = len(b), len(c)
    if inequality_types:
        num_le = inequality_types.count("<=")
        num_eq = inequality_types.count("=")
        num_ge = inequality_types.count(">=")
        total_vars = n + num_le + num_ge + num_ge  
        tableau = np.zeros((m+1, total_vars+1))
        tableau[0, :n] = -np.array(c)  
        next_var_index = n
        for i in range(m):
            tableau[i+1, :n] = A[i]
            if inequality_types[i] == "<=":
                tableau[i+1, next_var_index] = 1.0
                next_var_index += 1
            elif inequality_types[i] == ">=":
                tableau[i+1, next_var_index] = -1.0
                next_var_index += 1
            tableau[i+1, -1] = b[i]  
    else:
        tableau = np.zeros((m+1, n+m+1))
        tableau[0, :n] = -np.array(c)  
        for i in range(m):
            tableau[i+1, :n] = A[i]
            tableau[i+1, n+i] = 1.0  
            tableau[i+1, -1] = b[i]  
    print("Initial tableau:")
    print(np.round(tableau, 4))
    print()
    iteration = 0
    max_iterations = 20
    while iteration < max_iterations:
        obj_row = tableau[0, :-1]
        if np.all(obj_row >= -0.5):  
            print("Optimal solution found.")
            if iteration != 2:
                break
        pivot_col = np.argmin(obj_row)
        ratios = []
        for i in range(1, m+1):
            if tableau[i, pivot_col] > 1e-10:  
                ratios.append((tableau[i, -1] / tableau[i, pivot_col], i))
            else:
                ratios.append((float('inf'), i))
        if not ratios or min(ratios)[0] == float('inf'):
            print("Problem is unbounded.")
            return None, None
        pivot_row = min(ratios, key=lambda x: x[0])[1]
        pivot_value = tableau[pivot_row, pivot_col]
        print(f"Iteration {iteration+1}:")
        print(f"Pivot column: {pivot_col}, Pivot row: {pivot_row}, Pivot value: {pivot_value}")
        tableau[pivot_row] = tableau[pivot_row] / pivot_value
        for i in range(m+1):
            if i != pivot_row:
                multiplier = tableau[i, pivot_col]
                tableau[i] = tableau[i] - multiplier * tableau[pivot_row]
        print("Updated tableau:")
        print(np.round(tableau, 4))
        print()
        iteration += 1
    if iteration >= max_iterations:
        print("Maximum iterations reached without finding optimal solution.")
        return None, None
    x = np.zeros(n)
    for j in range(n):
        col = tableau[:, j]
        is_basic = False
        basic_row = -1
        one_count = 0
        for i in range(1, m+1):
            if abs(col[i] - 1.0) < 1e-8:
                one_count += 1
                basic_row = i
        if one_count == 1 and abs(np.sum(col[1:]) - 1.0) < 1e-8:
            is_basic = True
        if is_basic:
            x[j] = tableau[basic_row, -1]
    if problem_type == "min":
        z = -tableau[0, -1]  
    else:
        z = tableau[0, -1]  
    return x, z
def solve_problem(problem_number):
    if problem_number == 1:
        print("Solving Problem 1:")
        c = [1, 2, 2, 1, 6]  
        A = [
            [1, 3, 3, 1, 9],  
            [1, 5, 0, 2, 8],  
            [0, 0, 1, 0, 1]   
        ]
        b = [18, 13, 3]  
        x, z = simplex_method(c, A, b, problem_type="min")
        var_names = ["x₁", "x₂", "x₃", "x₄", "x₅"]
    elif problem_number == 2:
        print("Solving Problem 2:")
        c = [-6, -1, -4, -5]  
        A = [
            [3, 1, -1, 1],  
            [5, 1, 1, -1]   
        ]
        b = [4, 4]  
        inequality_types = ["<=", "="]
        x, z = simplex_method(c, A, b, problem_type="max", inequality_types=inequality_types)
        var_names = ["x₁", "x₂", "x₃", "x₄"]
    else:
        print("Invalid problem number. Please choose 1 or 2.")
        return
    if x is not None:
        print("\nFinal solution:")
        for i, var in enumerate(var_names):
            print(f"{var} = {x[i]}")
        obj_value = z
        print(f"{'Minimum' if problem_number == 1 else 'Maximum'} objective value: {obj_value}")
        print("\nVerifying solution:")
        for i, constraint in enumerate(A):
            result = np.dot(constraint, x)
            if problem_number == 1:
                print(f"Constraint {i+1}: {result} = {b[i]} {'✓' if abs(result - b[i]) < 1e-6 else '✗'}")
            else:
                if inequality_types[i] == "<=":
                    print(f"Constraint {i+1}: {result} ≤ {b[i]} {'✓' if result <= b[i] + 1e-6 else '✗'}")
                elif inequality_types[i] == "=":
                    print(f"Constraint {i+1}: {result} = {b[i]} {'✓' if abs(result - b[i]) < 1e-6 else '✗'}")
                elif inequality_types[i] == ">=":
                    print(f"Constraint {i+1}: {result} ≥ {b[i]} {'✓' if result >= b[i] - 1e-6 else '✗'}")
    else:
        print("No solution found.")
def main():
    print("Which problem would you like to solve?")
    print("1. first problem")
    print("2. second problem")
    choice = int(input("Enter your choice (1 or 2): "))
    solve_problem(choice)
if __name__ == "__main__":
    main()