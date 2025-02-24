import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog

inequalities = [
    (10, 3, 30, '>='),
    (-1, 1, 5, '<='),
    (1, 1, 10, '<='),
    (0, 1, 2, '>='),
]

# inequalities = [
#     (1, -1, 1, '<='),
#     (1, 1, 2, '<='),
#     (1, -2, 0, '>='),
#     (2, 2, 1, '>='),
# ]



A_ub = []
b_ub = []

for a1, a2, b, inequality_type in inequalities:
    if inequality_type == '<=':
        A_ub.append([a1, a2])
        b_ub.append(b)
    else:
        A_ub.append([-a1, -a2])
        b_ub.append(-b)

c = [1, 3]
# c = [1, -2]

res_min = linprog(c, A_ub=A_ub, b_ub=b_ub, method='highs')
res_max = linprog([-c[0], -c[1]], A_ub=A_ub, b_ub=b_ub, method='highs')

min_point = res_min.x if res_min.success else None
max_point = res_max.x if res_max.success else None
min_value = res_min.fun if res_min.success else None
max_value = -res_max.fun if res_max.success else None

x_values = np.linspace(-10, 10, 1000)
plt.figure(figsize=(8, 6))

lambdas = []
for a1, a2, b, inequality_type in inequalities:
    if a2 != 0:
        lambdas.append(lambda x2, a1=a1, a2=a2, b=b: (b - a1 * x2) / a2)
    else:
        lambdas.append(lambda x2, a1=a1, b=b: np.full_like(x2, b / a1))

for i, (a1, a2, b, inequality_type) in enumerate(inequalities):
    y_values = lambdas[i](x_values)
    
    if inequality_type == '<=':
        plt.plot(x_values, y_values, label=f'{a1}x₁ + {a2}x₂ <= {b}', color='blue')
        plt.fill_between(x_values, y_values, np.full_like(x_values, 10), color='blue', where=y_values <= 10, alpha=0.3)
    else:  # '>='
        plt.plot(x_values, y_values, label=f'{a1}x₁ + {a2}x₂ >= {b}', color='red')
        plt.fill_between(x_values, y_values, np.full_like(x_values, -10), color='red', where=y_values >= -10, alpha=0.3)



if min_point is not None:
    plt.scatter(min_point[0], min_point[1], color='red', s=100, 
                label=f'Мінімум: f={min_value:.2f}')
    plt.annotate(f'({min_point[0]:.2f}, {min_point[1]:.2f})', 
                (min_point[0], min_point[1]), xytext=(10, 10), 
                textcoords='offset points')
    
    if c[1] != 0:
        slope = -c[0]/c[1]
        x_min_line = np.array([min_point[0]-5, min_point[0]+5])
        y_min_line = slope * (x_min_line - min_point[0]) + min_point[1]
        plt.plot(x_min_line, y_min_line, 'r--', 
                label=f'Мін. лінія рівня: {c[0]}x₁ + {c[1]}x₂ = {min_value:.2f}')

if max_point is not None:
    plt.scatter(max_point[0], max_point[1], color='blue', s=100, 
                label=f'Максимум: f={max_value:.2f}')
    plt.annotate(f'({max_point[0]:.2f}, {max_point[1]:.2f})', 
                (max_point[0], max_point[1]), xytext=(10, 10), 
                textcoords='offset points')
    
    if c[1] != 0:
        slope = -c[0]/c[1]
        x_max_line = np.array([max_point[0]-5, max_point[0]+5])
        y_max_line = slope * (x_max_line - max_point[0]) + max_point[1]
        plt.plot(x_max_line, y_max_line, 'b--', 
                label=f'Макс. лінія рівня: {c[0]}x₁ + {c[1]}x₂ = {max_value:.2f}')

plt.xlabel('x₁')
plt.ylabel('x₂')
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.grid(True, linestyle='--', alpha=0.7)
plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
plt.legend()

plt.show()

print(f"Мінімум: f({min_point[0]:.2f}, {min_point[1]:.2f}) = {min_value:.2f}")
print(f"Максимум: f({max_point[0]:.2f}, {max_point[1]:.2f}) = {max_value:.2f}")