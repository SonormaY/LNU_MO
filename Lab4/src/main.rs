use std::io::{self, Write};
use std::collections::HashSet;

fn main() {
    // Constants
    let inf: i32 = 1000;
    
    println!("Transportation Problem Solver");
    println!("============================");
    
    // Ask user if they want to use default data or enter manually
    print!("Do you want to use default data? (y/n): ");
    io::stdout().flush().unwrap();
    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");
    let use_default = input.trim().to_lowercase() == "y";
    
    // Initialize variables
    let n: usize;
    let m: usize;
    let mut original_grid: Vec<Vec<i32>>;
    let supply: Vec<i32>;
    let demand: Vec<i32>;
    
    if use_default {
        // Use default data
        original_grid = vec![
            vec![4, 5, 6, 6],
            vec![6, 7, 4, 9],
            vec![7, 6, 8, 4]
        ];
        supply = vec![45, 15, 30];
        demand = vec![20, 22, 18, 30];
        
        n = original_grid.len();
        m = original_grid[0].len();
        
        println!("Using default data:");
        print_problem_state(&original_grid, &supply, &demand, 0);
    } else {
        // Get grid dimensions from user
        print!("Enter number of rows (supply sources): ");
        io::stdout().flush().unwrap();
        n = read_usize();
        
        print!("Enter number of columns (demand destinations): ");
        io::stdout().flush().unwrap();
        m = read_usize();
        
        // Input grid (cost matrix)
        println!("Enter the cost matrix ({} x {}):", n, m);
        original_grid = Vec::with_capacity(n);
        for i in 0..n {
            println!("Enter row {} (space-separated values):", i + 1);
            let row = read_vec_i32(m);
            original_grid.push(row);
        }
        
        // Input supply
        println!("Enter supply values ({} values):", n);
        supply = read_vec_i32(n);
        
        // Input demand
        println!("Enter demand values ({} values):", m);
        demand = read_vec_i32(m);
        
        // Print initial problem state
        println!("\nInitial Problem State:");
        print_problem_state(&original_grid, &supply, &demand, 0);
    }
    
    // Check if problem is balanced
    let total_supply: i32 = supply.iter().sum();
    let total_demand: i32 = demand.iter().sum();
    
    if total_supply != total_demand {
        println!("Warning: Problem is unbalanced. Total supply ({}) != Total demand ({}).", 
                 total_supply, total_demand);
        println!("Will attempt to solve anyway, but results may be suboptimal.");
    }
    
    // Ask user which method to use for initial solution
    print!("Choose initial solution method (1: Minimum Cost Method, 2: Vogel's Approximation Method): ");
    io::stdout().flush().unwrap();
    let mut method_input = String::new();
    io::stdin().read_line(&mut method_input).expect("Failed to read line");
    let method_choice = method_input.trim().parse::<usize>().unwrap_or(1);
    
    let mut allocations: Vec<(usize, usize, i32)>;
    let mut total_cost: i32;
    
    // Solve using the chosen method
    if method_choice == 2 {
        println!("\nSolving using Vogel's Approximation Method (VAM)...");
        (allocations, total_cost) = solve_vam(&original_grid, &supply, &demand);
    } else {
        println!("\nSolving using Minimum Cost Method...");
        (allocations, total_cost) = solve_minimum_cost(&original_grid, &supply, &demand);
    }
    
    // Print the initial solution
    println!("\nInitial Basic Feasible Solution:");
    print_allocations(&allocations, &original_grid, total_cost);
    
    // Ask if user wants to optimize using potential method
    print!("Do you want to optimize using Potential Method? (y/n): ");
    io::stdout().flush().unwrap();
    let mut opt_input = String::new();
    io::stdin().read_line(&mut opt_input).expect("Failed to read line");
    let optimize = opt_input.trim().to_lowercase() == "y";
    
    if optimize {
        println!("\nOptimizing solution using Potential Method...");
        let (optimized_allocations, optimized_cost) = optimize_potential_method(
            &original_grid, &allocations, n, m, &supply, &demand);
        
        // Print the optimized solution
        println!("\nOptimized Solution after Potential Method:");
        print_allocations(&optimized_allocations, &original_grid, optimized_cost);
        
        println!("\nImprovement: {} -> {} (Saved: {})", 
                 total_cost, optimized_cost, total_cost - optimized_cost);
    }
}

// Function to solve using Minimum Cost Method
fn solve_minimum_cost(grid: &Vec<Vec<i32>>, supply: &Vec<i32>, demand: &Vec<i32>) -> (Vec<(usize, usize, i32)>, i32) {
    let n = grid.len();
    let m = grid[0].len();
    let inf = 1000;
    
    // Create working copies
    let mut working_grid = grid.clone();
    let mut working_supply = supply.clone();
    let mut working_demand = demand.clone();
    
    // Initialize answer and allocations
    let mut ans = 0;
    let mut allocations: Vec<(usize, usize, i32)> = Vec::new(); // (row, col, amount)   
    let mut iteration = 1;
    
    // Loop runs until both demand and supply are exhausted
    while working_supply.iter().any(|&x| x > 0) && working_demand.iter().any(|&x| x > 0) {
        println!("\nIteration {}:", iteration);
        
        // Find minimum cost cell in entire grid
        let mut min_cost = inf;
        let mut min_i = 0;
        let mut min_j = 0;
        
        for i in 0..n {
            for j in 0..m {
                if working_grid[i][j] < min_cost && working_grid[i][j] != inf {
                    min_cost = working_grid[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }
        
        if min_cost == inf {
            println!("No more valid allocations possible. Problem might be unbalanced.");
            break;
        }
        
        println!("Selected cell with minimum cost {} at ({}, {})", 
                 min_cost, min_i + 1, min_j + 1);
        
        // Calculate minimum of supply and demand
        let allocation = std::cmp::min(working_supply[min_i], working_demand[min_j]);
        ans += allocation * min_cost;
        
        // Record this allocation
        allocations.push((min_i, min_j, allocation));
        
        // Subtract allocation from supply and demand
        working_supply[min_i] -= allocation;
        working_demand[min_j] -= allocation;
        
        println!("Allocated {} units from source {} to destination {} at cost {}", 
                 allocation, min_i + 1, min_j + 1, min_cost);
        
        // If demand is satisfied, eliminate the column
        if working_demand[min_j] == 0 {
            for r in 0..n {
                working_grid[r][min_j] = inf;
            }
            println!("Demand at destination {} is now satisfied", min_j + 1);
        }
        // If supply is exhausted, eliminate the row
        if working_supply[min_i] == 0 {
            for c in 0..m {
                working_grid[min_i][c] = inf;
            }
            println!("Supply at source {} is now exhausted", min_i + 1);
        }
        
        // Print state after this iteration
        print_problem_state(&working_grid, &working_supply, &working_demand, iteration);
        
        iteration += 1;
    }
    
    (allocations, ans)
}

// Function to solve using Vogel's Approximation Method (VAM)
fn solve_vam(grid: &Vec<Vec<i32>>, supply: &Vec<i32>, demand: &Vec<i32>) -> (Vec<(usize, usize, i32)>, i32) {
    let n = grid.len();
    let m = grid[0].len();
    let inf = 1000;
    
    // Create working copies
    let mut working_grid = grid.clone();
    let mut working_supply = supply.clone();
    let mut working_demand = demand.clone();
    
    // Initialize answer and allocations
    let mut ans = 0;
    let mut allocations: Vec<(usize, usize, i32)> = Vec::new(); // (row, col, amount)
    let mut iteration = 1;
    
    // Helper function for finding row and column differences
    fn find_diff(grid: &Vec<Vec<i32>>) -> (Vec<i32>, Vec<i32>) {
        let n = grid.len();
        let m = grid[0].len();
        
        // Calculate row differences
        let mut row_diff = Vec::with_capacity(n);
        for i in 0..n {
            let mut arr: Vec<i32> = grid[i].iter()
                .filter(|&&x| x != 1000)
                .cloned()
                .collect();
            
            if arr.len() >= 2 {
                arr.sort();
                row_diff.push(arr[1] - arr[0]);
            } else if arr.len() == 1 {
                row_diff.push(0);
            } else {
                row_diff.push(-1); // All elements are INF
            }
        }
        
        // Calculate column differences
        let mut col_diff = Vec::with_capacity(m);
        for col in 0..m {
            let mut arr: Vec<i32> = (0..n)
                .filter_map(|i| {
                    let val = grid[i][col];
                    if val != 1000 { Some(val) } else { None }
                })
                .collect();
            
            if arr.len() >= 2 {
                arr.sort();
                col_diff.push(arr[1] - arr[0]);
            } else if arr.len() == 1 {
                col_diff.push(0);
            } else {
                col_diff.push(-1); // All elements are INF
            }
        }
        
        (row_diff, col_diff)
    }
    
    // Loop runs until both demand and supply are exhausted
    while working_supply.iter().any(|&x| x > 0) && working_demand.iter().any(|&x| x > 0) {
        // Find row and column differences
        let (row, col) = find_diff(&working_grid);
        
        println!("\nIteration {}:", iteration);
        println!("Row differences: {:?}", row);
        println!("Column differences: {:?}", col);
        
        // Find maximum element in row difference array
        let max_row_diff = *row.iter().max().unwrap_or(&-1);
        
        // Find maximum element in column difference array
        let max_col_diff = *col.iter().max().unwrap_or(&-1);
        
        println!("Max row difference: {}", max_row_diff);
        println!("Max column difference: {}", max_col_diff);
        
        let mut allocation_made = false;
        
        // If row difference max element is greater than or equal to column difference max element
        if max_row_diff >= max_col_diff && max_row_diff >= 0 {
            for (ind, val) in row.iter().enumerate() {
                if *val == max_row_diff && working_supply[ind] > 0 {
                    // Find minimum element in grid row where maximum difference was found
                    let valid_elements: Vec<(usize, i32)> = working_grid[ind].iter()
                        .enumerate()
                        .filter(|(j, &x)| x != inf && working_demand[*j] > 0)
                        .map(|(i, &x)| (i, x))
                        .collect();
                    
                    if valid_elements.is_empty() {
                        continue;
                    }
                    
                    let (min_index, min_val) = valid_elements.iter()
                        .min_by_key(|(_, val)| val)
                        .unwrap();
                    
                    println!("Selected row {} with min cost {} at column {}", 
                             ind + 1, min_val, min_index + 1);
                    
                    // Calculate minimum of supply and demand
                    let allocation = std::cmp::min(working_supply[ind], working_demand[*min_index]);
                    ans += allocation * min_val;
                    
                    // Record this allocation
                    allocations.push((ind, *min_index, allocation));
                    
                    // Subtract allocation from supply and demand
                    working_supply[ind] -= allocation;
                    working_demand[*min_index] -= allocation;
                    
                    println!("Allocated {} units from source {} to destination {} at cost {}", 
                             allocation, ind + 1, min_index + 1, min_val);
                    
                    // If demand is satisfied, eliminate the column
                    if working_demand[*min_index] == 0 {
                        for r in 0..n {
                            working_grid[r][*min_index] = inf;
                        }
                        println!("Demand at destination {} is now satisfied", min_index + 1);
                    }
                    // If supply is exhausted, eliminate the row
                    if working_supply[ind] == 0 {
                        for c in 0..m {
                            working_grid[ind][c] = inf;
                        }
                        println!("Supply at source {} is now exhausted", ind + 1);
                    }
                    
                    allocation_made = true;
                    break;
                }
            }
        }
        // If column difference max element is greater
        if !allocation_made && max_col_diff >= 0 {
            for (ind, val) in col.iter().enumerate() {
                if *val == max_col_diff && working_demand[ind] > 0 {
                    // Find minimum element in grid column where maximum difference was found
                    let mut min_val = inf;
                    let mut min_index = 0;
                    let mut found = false;
                    
                    for j in 0..n {
                        if working_grid[j][ind] < min_val && working_grid[j][ind] != inf && working_supply[j] > 0 {
                            min_val = working_grid[j][ind];
                            min_index = j;
                            found = true;
                        }
                    }
                    
                    if !found {
                        continue;
                    }
                    
                    println!("Selected column {} with min cost {} at row {}", 
                             ind + 1, min_val, min_index + 1);
                    
                    // Calculate minimum of supply and demand
                    let allocation = std::cmp::min(working_supply[min_index], working_demand[ind]);
                    ans += allocation * min_val;
                    
                    // Record this allocation
                    allocations.push((min_index, ind, allocation));
                    
                    // Subtract allocation from supply and demand
                    working_supply[min_index] -= allocation;
                    working_demand[ind] -= allocation;
                    
                    println!("Allocated {} units from source {} to destination {} at cost {}", 
                             allocation, min_index + 1, ind + 1, min_val);
                    
                    // If demand is satisfied, eliminate the column
                    if working_demand[ind] == 0 {
                        for r in 0..n {
                            working_grid[r][ind] = inf;
                        }
                        println!("Demand at destination {} is now satisfied", ind + 1);
                    }
                    // If supply is exhausted, eliminate the row
                    if working_supply[min_index] == 0 {
                        for c in 0..m {
                            working_grid[min_index][c] = inf;
                        }
                        println!("Supply at source {} is now exhausted", min_index + 1);
                    }
                    
                    allocation_made = true;
                    break;
                }
            }
        }
        
        if !allocation_made {
            println!("No more valid allocations possible. Problem might be unbalanced.");
            break;
        }
        
        // Print state after this iteration
        print_problem_state(&working_grid, &working_supply, &working_demand, iteration);
        
        iteration += 1;
    }
    
    // Ensure we have enough allocations for the potential method
    ensure_non_degeneracy(&mut allocations, n, m);
    
    (allocations, ans)
}

// Function to ensure non-degeneracy of the solution
fn ensure_non_degeneracy(allocations: &mut Vec<(usize, usize, i32)>, n: usize, m: usize) {
    // For potential method to work properly, we need exactly n + m - 1 allocations
    let required_allocations = n + m - 1;
    
    // Check if we have enough allocations
    if allocations.len() < required_allocations {
        println!("\nWarning: Solution is degenerate. Adding zero allocations...");
        
        // Create a set of current allocated cells for quick lookup
        let mut allocated_cells = HashSet::new();
        for &(i, j, _) in allocations.iter() {
            allocated_cells.insert((i, j));
        }
        
        // Add zero allocations until we have n + m - 1 allocations
        for i in 0..n {
            for j in 0..m {
                if allocations.len() >= required_allocations {
                    break;
                }
                
                if !allocated_cells.contains(&(i, j)) {
                    allocations.push((i, j, 0));
                    allocated_cells.insert((i, j));
                    println!("Added zero allocation at ({}, {})", i + 1, j + 1);
                }
            }
            
            if allocations.len() >= required_allocations {
                break;
            }
        }
    }
}

// Function to optimize using potential method
fn optimize_potential_method(
    original_grid: &Vec<Vec<i32>>,
    initial_allocations: &Vec<(usize, usize, i32)>,
    n: usize,
    m: usize,
    supply: &Vec<i32>,
    demand: &Vec<i32>
) -> (Vec<(usize, usize, i32)>, i32) {
    let mut allocations = initial_allocations.clone();
    let mut iteration = 1;
    let mut total_cost = calculate_total_cost(original_grid, &allocations);
    
    loop {
        println!("\nPotential Method - Iteration {}:", iteration);
        
        // Step 1: Calculate potentials (u_i and v_j)
        let (u, v, success) = calculate_potentials(original_grid, &allocations, n, m);
        
        if !success {
            println!("Failed to calculate potentials. Solution might be degenerate.");
            break;
        }
        
        println!("Row potentials (u): {:?}", u);
        println!("Column potentials (v): {:?}", v);
        
        // Step 2: Calculate opportunity costs (delta_ij)
        let mut max_improvement = 0;
        let mut best_cell = (0, 0);
        
        for i in 0..n {
            for j in 0..m {
                if !is_allocated(&allocations, i, j) {
                    let c_ij = original_grid[i][j];
                    let delta_ij = c_ij - (u[i] + v[j]);
                    
                    // If delta_ij < 0, we can improve the solution
                    if delta_ij < 0 && delta_ij < max_improvement {
                        max_improvement = delta_ij;
                        best_cell = (i, j);
                    }
                }
            }
        }
        
        // If no improvement is possible, we have the optimal solution
        if max_improvement == 0 {
            println!("Current solution is optimal. No further improvement possible.");
            break;
        }
        
        println!("Best cell for improvement: ({}, {}) with delta = {}", 
                 best_cell.0 + 1, best_cell.1 + 1, max_improvement);
        
        // Step 3: Find the closed loop and determine the cells to update
        let loop_cells = find_loop(&allocations, best_cell, n, m);
        
        println!("Closed loop: {:?}", loop_cells.iter()
                 .map(|&(i, j)| (i + 1, j + 1))
                 .collect::<Vec<(usize, usize)>>());
        
        // Step 4: Find the minimum allocation amount in "minus" cells
        let mut min_allocation = std::i32::MAX;
        for (idx, &(i, j)) in loop_cells.iter().enumerate() {
            if idx % 2 == 1 { // Odd index means "minus" cell
                if let Some(&(_, _, amount)) = allocations.iter()
                    .find(|&&(r, c, _)| r == i && c == j) {
                    min_allocation = std::cmp::min(min_allocation, amount);
                }
            }
        }
        
        println!("Minimum allocation in loop: {}", min_allocation);
        
        // Step 5: Update allocations
        for (idx, &(i, j)) in loop_cells.iter().enumerate() {
            if idx == 0 { // Entering cell (new allocation)
                allocations.push((i, j, min_allocation));
                println!("Added allocation of {} units at ({}, {})", 
                         min_allocation, i + 1, j + 1);
            } else {
                // Find the allocation in our list
                for alloc in &mut allocations {
                    if alloc.0 == i && alloc.1 == j {
                        if idx % 2 == 0 { // Even index means "plus" cell
                            alloc.2 += min_allocation;
                            println!("Increased allocation at ({}, {}) by {} to {}", 
                                     i + 1, j + 1, min_allocation, alloc.2);
                        } else { // Odd index means "minus" cell
                            alloc.2 -= min_allocation;
                            println!("Decreased allocation at ({}, {}) by {} to {}", 
                                     i + 1, j + 1, min_allocation, alloc.2);
                            
                            // If allocation becomes zero, mark for removal
                            if alloc.2 == 0 {
                                *alloc = (usize::MAX, usize::MAX, 0);
                            }
                        }
                        break;
                    }
                }
            }
        }
        
        // Remove zero allocations
        allocations.retain(|&(i, j, _)| i != usize::MAX && j != usize::MAX);
        
        // Calculate new total cost
        let new_cost = calculate_total_cost(original_grid, &allocations);
        println!("New total cost: {} (Improved by: {})", new_cost, total_cost - new_cost);
        
        // Check if the solution is still feasible
        if !is_solution_feasible(&allocations, n, m, supply, demand) {
            println!("Warning: Solution is no longer feasible. Stopping optimization.");
            return (initial_allocations.clone(), total_cost);
        }
        
        total_cost = new_cost;
        iteration += 1;
        
        // Ensure we have the right number of allocations (non-degeneracy)
        ensure_non_degeneracy(&mut allocations, n, m);
        
        // Print current allocations
        print_allocations(&allocations, original_grid, total_cost);
        
        // Limit iterations to prevent infinite loops
        if iteration > 20 {
            println!("Maximum iterations reached. Stopping optimization.");
            break;
        }
    }
    
    (allocations, total_cost)
}

// Function to calculate u_i and v_j potentials
fn calculate_potentials(
    grid: &Vec<Vec<i32>>,
    allocations: &Vec<(usize, usize, i32)>,
    n: usize,
    m: usize
) -> (Vec<i32>, Vec<i32>, bool) {
    let mut u = vec![i32::MAX; n];
    let mut v = vec![i32::MAX; m];
    
    // Set u_0 = 0 (arbitrary)
    u[0] = 0;
    
    // Mark which potentials have been determined
    let mut u_determined = vec![false; n];
    let mut v_determined = vec![false; m];
    u_determined[0] = true;
    
    // Track progress
    let mut progress_made = true;
    let mut iterations = 0;
    let max_iterations = n * m; // Avoid infinite loops
    
    while progress_made && iterations < max_iterations {
        progress_made = false;
        iterations += 1;
        
        // For each allocation, try to determine potentials
        for &(i, j, _) in allocations {
            if u_determined[i] && !v_determined[j] {
                // If u_i is known, calculate v_j
                v[j] = grid[i][j] - u[i];
                v_determined[j] = true;
                progress_made = true;
            } else if !u_determined[i] && v_determined[j] {
                // If v_j is known, calculate u_i
                u[i] = grid[i][j] - v[j];
                u_determined[i] = true;
                progress_made = true;
            }
        }
    }
    
    // Check if all potentials were determined
    let all_determined = u_determined.iter().all(|&x| x) && v_determined.iter().all(|&x| x);
    
    (u, v, all_determined)
}

// Function to check if a cell is allocated
fn is_allocated(allocations: &Vec<(usize, usize, i32)>, i: usize, j: usize) -> bool {
    allocations.iter().any(|&(r, c, _)| r == i && c == j)
}

// Function to find a closed loop starting from the entering cell
fn find_loop(
    allocations: &Vec<(usize, usize, i32)>,
    entering_cell: (usize, usize),
    n: usize,
    m: usize
) -> Vec<(usize, usize)> {
    // First element is the entering cell
    let mut loop_cells = vec![entering_cell];
    
    // Create a grid to represent allocated cells
    let mut allocation_grid = vec![vec![false; m]; n];
    for &(i, j, _) in allocations {
        allocation_grid[i][j] = true;
    }
    
    // Recursively find the loop
    fn find_path(
        grid: &Vec<Vec<bool>>,
        start: (usize, usize),
        current: (usize, usize),
        path: &mut Vec<(usize, usize)>,
        visited_rows: &mut Vec<bool>,
        visited_cols: &mut Vec<bool>,
        n: usize,
        m: usize
    ) -> bool {
        // If we've returned to the starting point and the path has at least 4 cells, we found a loop
        if current == start && path.len() >= 4 {
            return true;
        }
        
        let (i, j) = current;
        
        // Try all rows in the current column j
        for row in 0..n {
            if row != i && grid[row][j] && !visited_rows[row] {
                visited_rows[row] = true;
                path.push((row, j));
                
                if find_path(grid, start, (row, j), path, visited_rows, visited_cols, n, m) {
                    return true;
                }
                
                visited_rows[row] = false;
                path.pop();
            }
        }
        
        // Try all columns in the current row i
        for col in 0..m {
            if col != j && grid[i][col] && !visited_cols[col] {
                visited_cols[col] = true;
                path.push((i, col));
                
                if find_path(grid, start, (i, col), path, visited_rows, visited_cols, n, m) {
                    return true;
                }
                
                visited_cols[col] = false;
                path.pop();
            }
        }
        
        // For the starting cell, we need to try empty rows/cols too
        if current == start {
            // Try potential horizontal moves
            for col in 0..m {
                if col != j {
                    path.push((i, col));
                    visited_cols[col] = true;
                    
                    if find_path(grid, start, (i, col), path, visited_rows, visited_cols, n, m) {
                        return true;
                    }
                    
                    visited_cols[col] = false;
                    path.pop();
                }
            }
            
            // Try potential vertical moves
            for row in 0..n {
                if row != i {
                    path.push((row, j));
                    visited_rows[row] = true;
                    
                    if find_path(grid, start, (row, j), path, visited_rows, visited_cols, n, m) {
                        return true;
                    }
                    
                    visited_rows[row] = false;
                    path.pop();
                }
            }
        }
        
        false
    }
    
    // Add entering cell to allocation grid temporarily
    allocation_grid[entering_cell.0][entering_cell.1] = true;
    
    // Find the loop
    let mut visited_rows = vec![false; n];
    let mut visited_cols = vec![false; m];
    let found = find_path(
        &allocation_grid,
        entering_cell,
        entering_cell,
        &mut loop_cells,
        &mut visited_rows,
        &mut visited_cols,
        n,
        m
    );
    
    if !found {
        // Simple fallback method: find a path through allocated cells
        loop_cells.clear();
        loop_cells.push(entering_cell);
        
        // Find row with allocation in the same column
        for &(i, j, _) in allocations {
            if j == entering_cell.1 {
                loop_cells.push((i, j));
                
                // Find column with allocation in the same row
                for &(r, c, _) in allocations {
                    if r == i && c != j {
                        loop_cells.push((r, c));
                        
                        // Find row with allocation in the same column
                        for &(p, q, _) in allocations {
                            if q == c && p == entering_cell.0 {
                                loop_cells.push((p, q));
                                return loop_cells;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Return the loop cells
    loop_cells
}

// Function to calculate total cost of a solution
fn calculate_total_cost(grid: &Vec<Vec<i32>>, allocations: &Vec<(usize, usize, i32)>) -> i32 {
    allocations.iter()
        .map(|&(i, j, amount)| grid[i][j] * amount)
        .sum()
}

// Function to check if the solution is feasible
fn is_solution_feasible(
    allocations: &Vec<(usize, usize, i32)>,
    n: usize,
    m: usize,
    supply: &Vec<i32>,
    demand: &Vec<i32>
) -> bool {
    // Check if all supply is satisfied
    let mut row_sum = vec![0; n];
    for &(i, _, amount) in allocations {
        row_sum[i] += amount;
    }
    
    for i in 0..n {
        if row_sum[i] != supply[i] {
            println!("Supply constraint violated at row {}: {} != {}", i + 1, row_sum[i], supply[i]);
            return false;
        }
    }
    
    // Check if all demand is satisfied
    let mut col_sum = vec![0; m];
    for &(_, j, amount) in allocations {
        col_sum[j] += amount;
    }
    
    for j in 0..m {
        if col_sum[j] != demand[j] {
            println!("Demand constraint violated at column {}: {} != {}", j + 1, col_sum[j], demand[j]);
            return false;
        }
    }
    
    true
}

// Function to print allocations
fn print_allocations(
    allocations: &Vec<(usize, usize, i32)>,
    grid: &Vec<Vec<i32>>,
    total_cost: i32
) {
    println!("\nAllocations:");
    println!("---------------------------------");
    println!("| Source | Dest | Units | Cost |");
    println!("---------------------------------");
    
    for &(row, col, amount) in allocations {
        if amount > 0 { // Only print non-zero allocations
            println!("| {:^6} | {:^4} | {:^5} | {:^4} |", 
                     row + 1, col + 1, amount, grid[row][col]);
        }
    }
    
    println!("---------------------------------");
    println!("Total transportation cost: {}", total_cost);
}

// Helper function to read a usize from stdin
fn read_usize() -> usize {
    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");
    input.trim().parse().expect("Please enter a valid number")
}

// Helper function to read a vector of i32 values
fn read_vec_i32(size: usize) -> Vec<i32> {
    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");
    
    let values: Vec<i32> = input
        .trim()
        .split_whitespace()
        .map(|s| s.parse().expect("Please enter valid numbers"))
        .collect();
    
    if values.len() != size {
        println!("Warning: Expected {} values, got {}. Using what was provided.", 
                 size, values.len());
    }
    
    values
}

// Function to print the current state of the problem
fn print_problem_state(grid: &Vec<Vec<i32>>, supply: &Vec<i32>, demand: &Vec<i32>, iteration: usize) {
    let n = grid.len();
    let m = grid[0].len();
    
    if iteration > 0 {
        println!("\nState after iteration {}:", iteration);
    }
    
    // Print cost matrix with supply values
    println!("Cost Matrix & Supply:");
    
    // Print header with destination numbers
    print!("     | ");
    for j in 0..m {
        print!("D{:^3} | ", j + 1);
    }
    println!("Supply");
    
    // Print separator
    print!("-----|");
    for _ in 0..m {
        print!("------|");
    }
    println!("-------");
    
    // Print each row with source number
    for i in 0..n {
        print!("S{:<3} | ", i + 1);
        for j in 0..m {
            if grid[i][j] == 1000 {
                print!("  X  | ");
            } else {
                print!("{:^4} | ", grid[i][j]);
            }
        }
        println!("{:^5}", supply[i]);
    }
    
    // Print demand row
    print!("Dmnd | ");
    for j in 0..m {
        print!("{:^4} | ", demand[j]);
    }
    println!();
    
    // Print separator
    print!("-----|");
    for _ in 0..m {
        print!("------|");
    }
    println!();
}
