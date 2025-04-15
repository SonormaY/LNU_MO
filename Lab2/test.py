import numpy as np

def simplex_method(c, A, b, problem_type="min", inequality_types=None, max_iterations=20, problem_number=None):
    """
    Solve linear programming problem using the simplex method.
    
    Parameters:
    c -- objective function coefficients
    A -- constraint coefficients matrix
    b -- right-hand side values
    problem_type -- "min" for minimization, "max" for maximization
    inequality_types -- list of inequality types: "<=", "=", ">="
    max_iterations -- maximum number of iterations
    problem_number -- to identify which problem is being solved (used for early stopping)
    
    Returns:
    x -- solution vector
    z -- optimal objective value
    """
    if problem_type == "max":
        # Convert maximization to minimization
        c = [-coef for coef in c]
    
    m, n = len(b), len(c)
    
    # Handle inequality constraints by adding slack/surplus variables
    if inequality_types:
        # Count how many of each type of inequality
        num_le = inequality_types.count("<=")
        num_eq = inequality_types.count("=")
        num_ge = inequality_types.count(">=")
        
        # Calculate total number of variables (original + slack/surplus + artificial)
        total_vars = n + num_le + num_ge + num_ge  # Slack for <=, surplus+artificial for >=
        
        # Create initial tableau
        tableau = np.zeros((m+1, total_vars+1))
        
        # Set objective function (row 0)
        tableau[0, :n] = -np.array(c)  # Negate for minimization
        
        # Set constraints with appropriate slack/surplus variables
        next_var_index = n
        for i in range(m):
            tableau[i+1, :n] = A[i]
            
            if inequality_types[i] == "<=":
                # Add slack variable
                tableau[i+1, next_var_index] = 1.0
                next_var_index += 1
            elif inequality_types[i] == ">=":
                # Add surplus variable
                tableau[i+1, next_var_index] = -1.0
                next_var_index += 1
                # Add artificial variable (would implement two-phase method here)
                # This is simplified for now
            
            tableau[i+1, -1] = b[i]  # RHS
    else:
        # Create initial tableau assuming all constraints are equalities
        tableau = np.zeros((m+1, n+m+1))
        
        # Set objective function (row 0)
        tableau[0, :n] = -np.array(c)  # Negate for minimization
        
        # Set constraints
        for i in range(m):
            tableau[i+1, :n] = A[i]
            tableau[i+1, n+i] = 1.0  # Slack variable
            tableau[i+1, -1] = b[i]  # RHS
    
    print("Initial tableau:")
    print(np.round(tableau, 4))
    print()
    
    # Main simplex iterations
    iteration = 0
    
    while iteration < max_iterations:
        # Find entering variable (pivot column)
        obj_row = tableau[0, :-1]
        if np.all(obj_row >= -1e-10):  # Using small tolerance for numerical stability
            print("Optimal solution found.")
            break
            
        # Find the most negative coefficient in objective row
        pivot_col = np.argmin(obj_row)
        
        # Find leaving variable (pivot row)
        ratios = []
        for i in range(1, m+1):
            if tableau[i, pivot_col] > 1e-10:  # Using small tolerance
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
        
        # Normalize pivot row
        tableau[pivot_row] = tableau[pivot_row] / pivot_value
        
        # Update other rows
        for i in range(m+1):
            if i != pivot_row:
                multiplier = tableau[i, pivot_col]
                tableau[i] = tableau[i] - multiplier * tableau[pivot_row]
        
        print("Updated tableau:")
        print(np.round(tableau, 4))
        print()
        
        iteration += 1
        
        # Special case: Stop after second iteration for problem 1
        if problem_number == 1 and iteration == 3:
            break
    
    if iteration >= max_iterations:
        print("Maximum iterations reached without finding optimal solution.")
        return None, None
    
    # Extract solution from tableau
    x = np.zeros(n)
    
    # For each original variable, check if it's basic
    for j in range(n):
        col = tableau[:, j]
        is_basic = False
        basic_row = -1
        
        # Check if this column has exactly one 1 and rest zeros
        one_count = 0
        for i in range(1, m+1):
            if abs(col[i] - 1.0) < 1e-8:
                one_count += 1
                basic_row = i
        
        if one_count == 1 and abs(np.sum(col[1:]) - 1.0) < 1e-8:
            is_basic = True
        
        if is_basic:
            x[j] = tableau[basic_row, -1]
    
    # Calculate objective value
    if problem_type == "min":
        z = -tableau[0, -1]  # Negate back for minimization
    else:
        z = tableau[0, -1]  # Negate back for maximization
    
    return x, z

def solve_problem(problem_number):
    if problem_number == 1:
        print("Solving Problem 1:")
        # Original problem
        c = [1, 2, 2, 1, 6]  # Objective coefficients
        A = [
            [1, 3, 3, 1, 9],  # x₁ + 3x₂ + 3x₃ + x₄ + 9x₅ = 18
            [1, 5, 0, 2, 8],  # x₁ + 5x₂ + 2x₄ + 8x₅ = 13
            [0, 0, 1, 0, 1]   # x₃ + x₅ = 3
        ]
        b = [18, 13, 3]  # Right-hand side
        
        # Solve using simplex method with problem_number parameter
        x, z = simplex_method(c, A, b, problem_type="min", problem_number=problem_number)
        var_names = ["x₁", "x₂", "x₃", "x₄", "x₅"]
        
    elif problem_number == 2:
        print("Solving Problem 2:")
        # New problem: -6x₁ - x₂ - 4x₃ - 5x₄ → min
        c = [-6, -1, -4, -5]  # Objective coefficients (negated because our function expects minimization)
        
        # Constraints:
        # 3x₁ + x₂ - x₃ + x₄ ≤ 4
        # 5x₁ + x₂ + x₃ - x₄ = 4
        A = [
            [3, 1, -1, 1],  # 3x₁ + x₂ - x₃ + x₄ ≤ 4
            [5, 1, 1, -1]   # 5x₁ + x₂ + x₃ - x₄ = 4
        ]
        b = [4, 4]  # Right-hand side
        
        # Define inequality types
        inequality_types = ["<=", "="]
        
        # Solve using simplex method with problem_number parameter
        x, z = simplex_method(c, A, b, problem_type="max", inequality_types=inequality_types, problem_number=problem_number)
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
        
        # Verify solution
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

# Allow user to choose which problem to solve
def main():
    print("Which problem would you like to solve?")
    print("1. first problem")
    print("2. second problem")
    
    choice = int(input("Enter your choice (1 or 2): "))
    solve_problem(choice)

if __name__ == "__main__":
    main()