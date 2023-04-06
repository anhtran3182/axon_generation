import numpy as np
import matplotlib.pyplot as plt
import generate_shape

from shapely.geometry import Point, Polygon, mapping

def generate_circle_center():
    x = np.random.uniform(min(boundary_x) + max_radius, max(boundary_x) - max_radius)
    y = np.random.uniform(min(boundary_y) + max_radius, max(boundary_y) - max_radius)    
    
    return Point(x, y)

def generate_circle_radius():
    return np.random.uniform(min_radius, max_radius)

def is_center_in_bound():
    # define the boundary and create a polygon
    polygon = Polygon(zip(boundary_x, boundary_x))

    center_point = generate_circle_center()

    while not polygon.contains(center_point):
        center_point = generate_circle_center()
        print(f"{center_point} is outside the polygon")
    
    print(f"{center_point} is outside the polygon")

def generate_circle_in_bound():
    polygon = Polygon(zip(boundary_x, boundary_x))

    min_x, min_y, max_x, max_y = polygon.bounds
    while True:
        point = Point(np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y))
        if polygon.contains(point):
            return np.array([point.x, point.y])

def circle_overlap(c1, c2):
    d = np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)
    return d < c1[2] + c2[2]
    
# define the shape's boundary as a list of x and y points
rad = 0.2
edgy = 0.05
a = generate_shape.get_random_points(n=7, scale=1) 
boundary_x, boundary_y, _ = generate_shape.get_bezier_curve(a,rad=rad, edgy=edgy)

# define requirements
n_circles = 10
min_radius = 1
max_radius = 3

print(is_center_in_bound())
