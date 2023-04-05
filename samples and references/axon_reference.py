'''
This file: generate axons from shapes that are generated randomly
Reference: https://www.matecdev.com/posts/random-points-in-polygon.html
'''

import geopandas as gpd
import numpy as np
import generate_shape
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, mapping

def generate_circle_centers():
    minx, miny, maxx, maxy = polygon.bounds
    while True:
        point = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
        if polygon.contains(point):
            return point

def generate_circle_radius():
    return np.random.uniform(min_radius, max_radius)

def generate_circles():
    circles = []

    while len(circles) < n_circles:
        center = generate_circle_centers()
        radius = generate_circle_radius()
        circle = center.buffer(radius)
        if not any(existing_circle.intersects(circle) for existing_circle in circles) and polygon.contains(circle):
            circles.append(circle)

    return circles

############# Boundary Set Up #############
# define the shape's boundary as a list of x and y points
rad = 0.2
edgy = 0.05
a = generate_shape.get_random_points(n=7, scale=1) 
boundary_x, boundary_y, _ = generate_shape.get_bezier_curve(a,rad=rad, edgy=edgy)

# make the polygon based on list
coords = list(zip(boundary_x, boundary_y))
polygon = Polygon(coords)


print(boundary_x)
print(boundary_y)

# Plot the polygon
xp,yp = polygon.exterior.xy
plt.plot(xp,yp)
    #plt.show()

# generate random points
points = generate_circle_centers(polygon, 20)

# define requirements
n_circles = 50
min_radius = 0.005
max_radius = 0.01
'''
# redefine the polygon allowing the generated circles to not touch the boundary 
# and are completely within the polygon
buffer_distance = min_radius * 0.5
polygon = polygon.buffer(-buffer_distance)

xp,yp = polygon.exterior.xy
plt.plot(xp,yp)

circles = generate_circles()
for circle in circles:
    xp,yp = circle.exterior.xy
    plt.plot(xp,yp)

plt.show()
'''

