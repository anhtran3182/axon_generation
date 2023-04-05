
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

import generate_shape

def generate_non_overlapping_circles(x, y, n_circles, min_radius, max_radius):
    # Create a path object that represents the shape
    path = create_path(x, y)

    # Generate circles that are contained within the shape and do not overlap
    circles = []
    while len(circles) < n_circles:
        # Generate a random circle center
        center = generate_random_center(x, y)

        # Check if the circle center is within the shape
        if path.contains_point(center):
            radius = generate_random_radius(min_radius, max_radius)
            circle = (center[0], center[1], radius)

            # Check if the circle overlaps with any existing circle
            if not any(circle_overlap(circle, c) for c in circles):
                circles.append(circle)
    return circles

def create_path(x, y):
    vertices = np.column_stack((x, y))
    codes = [Path.LINETO] * (len(vertices) - 1)
    codes[0] = Path.MOVETO
    return Path(vertices, codes)

def generate_random_center(x, y):
    return np.random.uniform(np.min(x), np.max(x), size=2)

def generate_random_radius(min_radius, max_radius):
    return np.random.uniform(min_radius, max_radius)

def circle_overlap(c1, c2):
    d = np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)
    return d < c1[2] + c2[2]

# Define the x and y coordinates of the arbitrary shape
#x = [1, 4, 4, 1, 0]
#y = [1, 1, 4, 4, 2]

rad = 0.2
edgy = 0.05

a = generate_shape.get_random_points(n=7, scale=1) 
x,y, _ = generate_shape.get_bezier_curve(a,rad=rad, edgy=edgy)

# Generate the circles within the arbitrary shape
n_circles = 20
min_radius = 0.5
max_radius = 2.5
circles = generate_non_overlapping_circles(x, y, n_circles, min_radius, max_radius)

# Plot the arbitrary shape and the circles
fig, ax = plt.subplots()
ax.plot(x, y)

for circle in circles:
    circle2 = plt.Circle(circle[:2], circle[2], color='blue', alpha=0.5)
    ax.add_artist(circle2)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
