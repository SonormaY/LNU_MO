use std::io::{self, Write};

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
    let mut grid: Vec<Vec<i32>>;
    let supply: Vec<i32>;
    let demand: Vec<i32>;
    
    if use_default {
        // Use default data
        grid = vec![
            vec![4, 5, 6, 6],
            vec![6, 7, 4, 9],
            vec![7, 6, 8, 4]
        ];
        supply = vec![45, 15, 30];
        demand = vec![20, 22, 18, 30];
        
        n = grid.len();
        m = grid[0].len();
        
        println!("Using default data:");
        print_problem_state(&grid, &supply, &demand, 0);
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
        grid = Vec::with_capacity(n);
        for i in 0..n {
            println!("Enter row {} (space-separated values):", i + 1);
            let row = read_vec_i32(m);
            grid.push(row);
        }
        
        // Input supply
        println!("Enter supply values ({} values):", n);
        supply = read_vec_i32(n);
        
        // Input demand
        println!("Enter demand values ({} values):", m);
        demand = read_vec_i32(m);
        
        // Print initial problem state
        println!("\nInitial Problem State:");
        print_problem_state(&grid, &supply, &demand, 0);
    }
    
    // Create working copies
    let mut working_grid = grid.clone();
    let mut working_supply = supply.clone();
    let mut working_demand = demand.clone();
    
    // Initialize answer
    let mut ans = 0;
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
    
    // Track allocations for pretty printing
    let mut allocations: Vec<(usize, usize, i32, i32)> = Vec::new(); // (row, col, amount, cost)
    
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
                if *val == max_row_diff {
                    // Find minimum element in grid row where maximum difference was found
                    let valid_elements: Vec<(usize, i32)> = working_grid[ind].iter()
                        .enumerate()
                        .filter(|(_, &x)| x != inf)
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
                    allocations.push((ind, *min_index, allocation, *min_val));
                    
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
                    else {
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
        else if max_col_diff >= 0 {
            for (ind, val) in col.iter().enumerate() {
                if *val == max_col_diff {
                    // Find minimum element in grid column where maximum difference was found
                    let mut min_val = inf;
                    let mut min_index = 0;
                    for j in 0..n {
                        if working_grid[j][ind] < min_val && working_grid[j][ind] != inf {
                            min_val = working_grid[j][ind];
                            min_index = j;
                        }
                    }
                    
                    if min_val == inf {
                        continue;
                    }
                    
                    println!("Selected column {} with min cost {} at row {}", 
                             ind + 1, min_val, min_index + 1);
                    
                    // Calculate minimum of supply and demand
                    let allocation = std::cmp::min(working_supply[min_index], working_demand[ind]);
                    ans += allocation * min_val;
                    
                    // Record this allocation
                    allocations.push((min_index, ind, allocation, min_val));
                    
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
                    else {
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
        
        // Print current allocations
        println!("\nCurrent Allocations:");
        println!("---------------------------------");
        println!("| Source | Dest | Units | Cost |");
        println!("---------------------------------");
        for (row, col, amount, cost) in &allocations {
            println!("| {:^6} | {:^4} | {:^5} | {:^4} |", 
                     row + 1, col + 1, amount, cost);
        }
        println!("---------------------------------");
        println!("Current total cost: {}", ans);
        
        iteration += 1;
    }
    
    // Print final solution
    println!("\n===============================");
    println!("Final Basic Feasible Solution");
    println!("===============================");
    println!("Allocations:");
    println!("---------------------------------");
    println!("| Source | Dest | Units | Cost |");
    println!("---------------------------------");
    for (row, col, amount, cost) in &allocations {
        println!("| {:^6} | {:^4} | {:^5} | {:^4} |", 
                 row + 1, col + 1, amount, cost);
    }
    println!("---------------------------------");
    println!("Total transportation cost: {}", ans);
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
