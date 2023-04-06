import numpy as np
import matplotlib.pyplot as plt

import generate_shape

def generate_circle_within_boundary(boundary_x, boundary_y, min_radius, max_radius, circles):
    # generate a single circle within the boundary
    while True:
        # generate random x and y coordinates within the boundary
        x = np.random.uniform(min(boundary_x) + max_radius, max(boundary_x) - max_radius)
        y = np.random.uniform(min(boundary_y) + max_radius, max(boundary_y) - max_radius)
        
        # generate random radius
        radius = np.random.uniform(min_radius, max_radius)
        
        # check if the circle is overbound
        if min(abs(boundary_x - x)) < radius or min(abs(boundary_y - y)) < radius:
            continue
        
        # check if the circle overlaps with any of the previously generated circles
        if not check_circle_overlaps(x, y, radius, circles):
            # return the circle
            return (x, y, radius)

def check_circle_overlaps(x, y, radius, circles):
    # check if a circle with center (x, y) and radius `radius` overlaps with any of the circles in `circles`
    for c in circles:
        if np.sqrt((x - c[0])**2 + (y - c[1])**2) < radius + c[2]:
            return True
    return False

def generate_circles_within_boundary(boundary_x, boundary_y, min_radius, max_radius, num_circles):
    # generate multiple circles within the boundary
    circles = []
    while len(circles) < num_circles:
        circle = generate_circle_within_boundary(boundary_x, boundary_y, min_radius, max_radius, circles)
        circles.append(circle)
    return circles

def plot_circles_within_boundary(boundary_x, boundary_y, circles):
    # plot the circles within the boundary
    fig, ax = plt.subplots()
    for c in circles:
        circle = plt.Circle((c[0], c[1]), c[2], fill=False)
        ax.add_artist(circle)
    ax.plot(boundary_x, boundary_y)
    plt.axis('equal')
    plt.show()

# define the shape's boundary as a list of x and y points
rad = 0.2
edgy = 0.05

a = generate_shape.get_random_points(n=7, scale=1) 
boundary_x, boundary_y, _ = generate_shape.get_bezier_curve(a,rad=rad, edgy=edgy)

# define the minimum and maximum radius of the circles
min_radius = 1
max_radius = 10

# generate circles
num_circles = 10
circles = generate_circles_within_boundary(boundary_x, boundary_y, min_radius, max_radius, num_circles)

# plot the circles
plot_circles_within_boundary(boundary_x, boundary_y, circles)
