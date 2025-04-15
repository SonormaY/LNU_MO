use plotters::prelude::*;
use std::f64;

// Define our function f(x) = 2x^2 - x + 2
fn f(x: f64) -> f64 {
    2.0 * x.powi(2) - x + 2.0
}

// Derivative of f(x)
fn df(x: f64) -> f64 {
    4.0 * x - 1.0
}

// DSC method to find initial interval [a,b]
fn dsc(x0: f64, h: f64) -> [f64; 2] {
    let mut x0 = x0;
    let mut h = h;
    let mut x1 = x0 + h;
    let mut f0 = f(x0);
    let mut f1 = f(x1);
    
    println!("DSC initialization:");
    println!("x0 = {}, f(x0) = {}", x0, f0);
    println!("x1 = {}, f(x1) = {}", x1, f1);
    
    // If f(x1) > f(x0), change direction
    if f1 > f0 {
        h = -h;
        x1 = x0 + h;
        f1 = f(x1);
        println!("Changed direction: h = {}", h);
        println!("New x1 = {}, f(x1) = {}", x1, f1);
    }
    
    println!("\nDSC iterations:");
    // While we're descending
    let mut iteration = 1;
    while f1 < f0 {
        println!("Iteration {}: x0 = {}, f(x0) = {}, x1 = {}, f(x1) = {}", 
                 iteration, x0, f0, x1, f1);
        h = h * 2.0;
        let x2 = x1 + h;
        let f2 = f(x2);
        
        println!("  h = {}, x2 = {}, f(x2) = {}", h, x2, f2);
        
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
        
        iteration += 1;
    }
    
    // Return interval [a, b] where a < b
    let mut interval = [x0, x1];
    if interval[0] > interval[1] {
        interval.swap(0, 1);
    }
    
    println!("\nFinal interval [a,b] = [{}, {}]", interval[0], interval[1]);
    interval
}

// Find the maximum of the absolute values of derivatives at critical and boundary points
fn find_lipschitz_constant(a: f64, b: f64) -> f64 {
    // Solve df(x) = 0 for critical points
    let critical_point = 1.0 / 4.0; // analytically derived from 4x - 1 = 0
    
    let boundary_points = [a, b];
    let mut all_points = vec![critical_point];
    all_points.extend_from_slice(&boundary_points);
    
    // Calculate derivative at each point
    let derivatives: Vec<f64> = all_points.iter().map(|&x| df(x).abs()).collect();
    
    println!("Critical point: x = {}", critical_point);
    println!("Derivatives at points: {:?}", derivatives);
    
    // Return maximum absolute derivative value
    *derivatives.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
}

// Get the value of the lower approximating function r_i(x) at point x
fn lower_approximation(x: f64, xk: f64, l: f64) -> f64 {
    f(xk) - l * (x - xk).abs()
}

// Find the maximum of the lower approximation function
fn find_max_of_lower_approx(a: f64, b: f64, points_with_values: &[(f64, f64)], l: f64) -> f64 {
    let num_samples = 1000;
    let step = (b - a) / (num_samples as f64);
    
    let mut max_x = a;
    let mut max_val = f64::NEG_INFINITY;
    
    for i in 0..=num_samples {
        let x = a + i as f64 * step;
        
        // Calculate r_i(x) - take maximum of all previous lower approximations
        let mut r_val = f64::NEG_INFINITY;
        for &(xk, fk) in points_with_values {
            let approx = fk - l * (x - xk).abs();
            if approx > r_val {
                r_val = approx;
            }
        }
        
        if r_val > max_val {
            max_val = r_val;
            max_x = x;
        }
    }
    
    max_x
}

// Plot function and lower approximations
fn plot_iteration(iteration: usize, x_values: &[f64], f_values: &[f64], 
                  r_values: &[f64], xk: f64, a: f64, b: f64) -> Result<(), Box<dyn std::error::Error>> {
    let filename = format!("broken_line_iteration_{}.png", iteration);
    let root = BitMapBackend::new(&filename, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let min_y = f_values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_y = f_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let margin = (max_y - min_y) * 0.1;
    
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Broken Line Method - Iteration {}", iteration), ("sans-serif", 30).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(a-0.5..b+0.5, min_y-margin..max_y+margin)?;
    
    chart.configure_mesh()
        .x_labels(10)
        .y_labels(10)
        .x_desc("x")
        .y_desc("y")
        .draw()?;
    
    // Draw the original function
    chart.draw_series(LineSeries::new(
        x_values.iter().zip(f_values.iter()).map(|(&x, &y)| (x, y)),
        &RED,
    ))?
    .label("f(x) = 2x^2 - x + 2")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    // Draw the lower approximation
    chart.draw_series(LineSeries::new(
        x_values.iter().zip(r_values.iter()).map(|(&x, &y)| (x, y)),
        &BLUE,
    ))?
    .label("Lower approximation r_i(x)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    
    // Mark the current point xk
    chart.draw_series(PointSeries::of_element(
        vec![(xk, f(xk))],
        5,
        &GREEN,
        &|c, s, st| {
            return EmptyElement::at(c)
                + Circle::new((0, 0), s, st.filled())
                + Text::new(format!("xk = {:.4}", xk), (10, 0), ("sans-serif", 15).into_font());
        },
    ))?;
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    println!("Plot saved as {}", filename);
    Ok(())
}

// Main broken line method algorithm
fn broken_line_method(x0: f64, h: f64, epsilon: f64) -> f64 {
    println!("Starting Broken Line Method with:");
    println!("x0 = {}, h = {}, epsilon = {}", x0, h, epsilon);
    
    // Step 1: Use DSC method to find interval [a,b]
    let [a, b] = dsc(x0, h);
    
    // Step 2: Find Lipschitz constant L
    let l = find_lipschitz_constant(a, b);
    println!("\nLipschitz constant L = {}", l);
    
    // Step 3: Initialize the broken line method
    let mut xk = (a + b) / 2.0; // Start from middle point
    let mut points_with_values = vec![(xk, f(xk))];
    let mut iteration = 1;
    
    println!("\nBroken Line Method iterations:");
    println!("Initial point: xk = {}, f(xk) = {}", xk, f(xk));
    
    // For visualization
    let num_plot_points = 500;
    let plot_x_values: Vec<f64> = (0..num_plot_points)
        .map(|i| a - 0.5 + (b - a + 1.0) * (i as f64) / (num_plot_points as f64 - 1.0))
        .collect();
    let plot_f_values: Vec<f64> = plot_x_values.iter().map(|&x| f(x)).collect();
    
    // Continue until the accuracy condition is met
    loop {
        // Calculate current lower approximation values for plotting
        let plot_r_values: Vec<f64> = plot_x_values.iter().map(|&x| {
            points_with_values.iter()
                .map(|&(xi, fi)| fi - l * (x - xi).abs())
                .fold(f64::NEG_INFINITY, f64::max)
        }).collect();
        
        // Try to plot this iteration
        if let Err(e) = plot_iteration(
            iteration, 
            &plot_x_values, 
            &plot_f_values, 
            &plot_r_values, 
            xk, 
            a, 
            b
        ) {
            println!("Error plotting iteration {}: {}", iteration, e);
        }
        
        // Find the maximum of the lower approximation to get the next point
        let next_xk = find_max_of_lower_approx(a, b, &points_with_values, l);
        
        // Calculate function value and lower approximation value at next_xk
        let f_next = f(next_xk);
        let r_next = points_with_values.iter()
            .map(|&(xi, fi)| fi - l * (next_xk - xi).abs())
            .fold(f64::NEG_INFINITY, f64::max);
        
        println!("Iteration {}: xk = {:.6}, f(xk) = {:.6}, r(xk) = {:.6}", 
                 iteration, next_xk, f_next, r_next);
        
        // Check the termination condition
        if f_next - r_next <= epsilon {
            println!("\nConvergence achieved!");
            println!("Final solution: x* = {:.6}, f(x*) = {:.6}", next_xk, f_next);
            return next_xk;
        }
        
        // Add the new point to our collection
        points_with_values.push((next_xk, f_next));
        xk = next_xk;
        iteration += 1;
        
        // Emergency break to prevent infinite loops
        if iteration > 20 {
            println!("Maximum iterations reached without convergence.");
            break;
        }
    }
    
    xk
}

fn main() {
    let x0 = 3.0;
    let h = 2.0;
    let epsilon = 0.01;
    
    println!("Finding minimum of f(x) = 2x^2 - x + 2");
    println!("Using the Broken Line Method (метод ламаних)");
    println!("------------------------------------------------");
    
    // Use the analytical solution as reference
    let analytical_min_x = 0.25; // solving 4x - 1 = 0
    let analytical_min_f = f(analytical_min_x);
    println!("Analytical minimum: x* = {}, f(x*) = {}", analytical_min_x, analytical_min_f);
    println!("------------------------------------------------\n");
    
    // Run the broken line method
    let result = broken_line_method(x0, h, epsilon);
    
    println!("\n------------------------------------------------");
    println!("Final result:");
    println!("Minimum found at x* = {:.6}", result);
    println!("Function value f(x*) = {:.6}", f(result));
    println!("True minimum at x = {:.6}", analytical_min_x);
    println!("Error: |x* - true_x*| = {:.6}", (result - analytical_min_x).abs());
}
