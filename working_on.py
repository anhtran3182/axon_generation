import numpy as np
import matplotlib.pyplot as plt
import generate_shape

from shapely.geometry import Point, Polygon

def generate_non_overlapping_circles(x, y, n_circles, min_radius, max_radius):
    # Create a path object that represents the shape
    coords = list(zip(boundary_x, boundary_x))
    polygon = Polygon(coords)

    point = generate_point(boundary_x, boundary_y, max_radius)
    
    # Generate circles that are contained within the shape and do not overlap
    circles = []
    while len(circles) < n_circles:
        # Check if the circle center is within the shape
        if polygon.contains(point):
            radius = np.random.uniform(min_radius, max_radius)
            circle = (center[0], center[1], radius)
            print(f"{point} is inside the polygon")     
            # Check if the circle overlaps with any existing circle
            if not any(circle_overlap(circle, c) for c in circles):
                circles.append(circle)
        
    return circles

def generate_point(boundary_x, boundary_y, max_radius):
    return np.random.uniform(np.min(boundary_x), np.max(boundary_x), size=2)

def circle_overlap(c1, c2):
    d = np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)
    return d < c1[2] + c2[2]
    
# define the shape's boundary as a list of x and y points
rad = 0.2
edgy = 0.05

a = generate_shape.get_random_points(n=7, scale=1) 
boundary_x, boundary_y, _ = generate_shape.get_bezier_curve(a,rad=rad, edgy=edgy)


# Generate the circles
center = np.array([0, 0])
radius = 5.0
n_circles = 10
min_radius = 0.5
max_radius = 1

generate_non_overlapping_circles(boundary_x, boundary_y, n_circles, min_radius, max_radius)