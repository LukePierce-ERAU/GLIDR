import numpy as np
import matplotlib.pyplot as plt

def load_airfoil_points(filename):
    points = []
    with open(filename, 'r') as f:
        for line in f:
            try:
                x, y = map(float, line.strip().split())
                points.append((x, y))
            except ValueError:
                continue  # skip headers or bad lines
    return np.array(points)

def compute_curve_length(points):
    diffs = np.diff(points, axis=0)
    segment_lengths = np.linalg.norm(diffs, axis=1)
    return np.sum(segment_lengths)

def find_stagnation_point(points):
    # Heuristic: find point closest to leading edge (minimum x)
    # then select the one closest to (0,0) from those
    min_x = np.min(points[:, 0])
    near_leading_edge = points[np.isclose(points[:, 0], min_x, atol=1e-3)]
    if len(near_leading_edge) == 0:
        return points[np.argmin(points[:, 0])]  # fallback
    distances = np.linalg.norm(near_leading_edge, axis=1)
    return near_leading_edge[np.argmin(distances)]

def plot_airfoil(points, stagnation_point=None):
    plt.figure(figsize=(10, 4))
    plt.plot(points[:, 0], points[:, 1], '-k', label='Airfoil')
    if stagnation_point is not None:
        plt.plot(stagnation_point[0], stagnation_point[1], 'ro', label='Stagnation Point')
        plt.text(stagnation_point[0], stagnation_point[1], '  stagnation', color='red')
    plt.axis('equal')
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Airfoil Shape with Stagnation Point')
    plt.legend()
    plt.tight_layout()
    plt.show()

# === MAIN SCRIPT ===
filename = 'airfoil.dat'
points = load_airfoil_points(filename)
length = compute_curve_length(points)
stagnation_point = find_stagnation_point(points)

print(f"Total curve length: {length:.6f}")
print(f"Stagnation point: {stagnation_point}")

plot_airfoil(points, stagnation_point)

