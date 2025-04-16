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

// Plot function and tangent line
fn plot_iteration(iteration: usize, x0: f64, x1: f64, a: f64, b: f64) -> Result<(), Box<dyn std::error::Error>> {
    let filename = format!("newton_iteration_{}.png", iteration);
    let root = BitMapBackend::new(&filename, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Calculate function values for plotting
    let num_points = 500;
    let step = (b - a) / (num_points as f64);
    let mut x_values = Vec::with_capacity(num_points);
    let mut f_values = Vec::with_capacity(num_points);
    
    for i in 0..num_points {
        let x = a + (i as f64) * step;
        x_values.push(x);
        f_values.push(f(x));
    }
    
    let min_y = *f_values.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let max_y = *f_values.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let margin = (max_y - min_y) * 0.1;
    
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Newton's Method - Iteration {}", iteration), ("sans-serif", 30).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(a..b, min_y-margin..max_y+margin)?;
    
    chart.configure_mesh()
        .x_labels(10)
        .y_labels(10)
        .x_desc("x")
        .y_desc("y")
        .draw()?;
    
    // Draw the original function
    chart.draw_series(LineSeries::new(
        x_values.iter().zip(f_values.iter()).map(|(&x, &y)| (x, y)),
        &BLUE,
    ))?
    .label("f(x) = 2x^2 - x + 2")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    
    // Draw the tangent line
    let f_x0 = f(x0);
    let df_x0 = df(x0);
    
    // Generate points for tangent line
    let tangent_x_min = a;
    let tangent_x_max = b;
    let tangent_y_min = f_x0 + df_x0 * (tangent_x_min - x0);
    let tangent_y_max = f_x0 + df_x0 * (tangent_x_max - x0);
    
    chart.draw_series(LineSeries::new(
        vec![(tangent_x_min, tangent_y_min), (tangent_x_max, tangent_y_max)],
        &RED,
    ))?
    .label("Tangent line")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    // Mark current point x0
    chart.draw_series(PointSeries::of_element(
        vec![(x0, f_x0)],
        5,
        &GREEN,
        &|c, s, st| {
            return EmptyElement::at(c)
                + Circle::new((0, 0), s, st.filled())
                + Text::new(format!("x0 = {:.4}", x0), (10, 0), ("sans-serif", 15).into_font());
        },
    ))?;
    
    // Mark next point x1
    chart.draw_series(PointSeries::of_element(
        vec![(x1, f(x1))],
        5,
        &MAGENTA,
        &|c, s, st| {
            return EmptyElement::at(c)
                + Circle::new((0, 0), s, st.filled())
                + Text::new(format!("x1 = {:.4}", x1), (10, 0), ("sans-serif", 15).into_font());
        },
    ))?;
    
    // Draw line from x1 to function
    chart.draw_series(LineSeries::new(
        vec![(x1, 0.0), (x1, f(x1))],
        &BLACK.mix(0.5),
    ))?;
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    println!("Plot saved as {}", filename);
    Ok(())
}

// Newton's method implementation
fn newton_method(x0: f64, epsilon: f64, max_iterations: usize) -> f64 {
    let mut x0 = x0;
    let plot_range = [-1.0, 1.0]; // Plot range [a, b]
    
    println!("Starting Newton's Method with:");
    println!("x0 = {}, epsilon = {}", x0, epsilon);
    println!("\nNewton's Method iterations:");
    
    for i in 0..max_iterations {
        let fx0 = f(x0);
        let fprime_x0 = df(x0);
        
        println!("Iteration {}: x = {:.6}, f(x) = {:.6}, f'(x) = {:.6}", 
                 i+1, x0, fx0, fprime_x0);
        
        // Newton's formula for the next approximation
        let x1 = x0 - fx0 / fprime_x0;
        
        // Plot this iteration
        if let Err(e) = plot_iteration(i+1, x0, x1, plot_range[0], plot_range[1]) {
            println!("Error plotting iteration {}: {}", i+1, e);
        }
        
        // Check convergence
        if (x1 - x0).abs() < epsilon {
            println!("\nConvergence achieved after {} iterations!", i+1);
            println!("Final solution: x* = {:.6}, f(x*) = {:.6}", x1, f(x1));
            return x1;
        }
        
        // Update x0 for next iteration
        x0 = x1;
    }
    
    println!("Maximum iterations reached without convergence.");
    x0
}

fn main() {
    let x0 = 1.0;
    let epsilon = 0.01;
    let max_iterations = 5;
    
    println!("Finding minimum of f(x) = 2x^2 - x + 2");
    println!("Using Newton's Method (метод дотичних)");
    println!("------------------------------------------------");
    
    // Use the analytical solution as reference
    let analytical_min_x = 0.25; // solving 4x - 1 = 0
    let analytical_min_f = f(analytical_min_x);
    println!("Analytical minimum: x* = {}, f(x*) = {}", analytical_min_x, analytical_min_f);
    println!("------------------------------------------------\n");
    
    // Run Newton's method
    let result = newton_method(x0, epsilon, max_iterations);
    
    println!("\n------------------------------------------------");
    println!("Final result:");
    println!("Minimum found at x* = {:.6}", result);
    println!("Function value f(x*) = {:.6}", f(result));
    println!("True minimum at x = {:.6}", analytical_min_x);
    println!("Error: |x* - true_x*| = {:.6}", (result - analytical_min_x).abs());
}
