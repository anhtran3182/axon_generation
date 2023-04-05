import geopandas as gpd
import numpy as np
import generate_shape
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, mapping

def generate_circle_centers(number):
    points = []
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < number:
        pnt = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
    
    return points

def generate_circle_radius():
    return np.random.uniform(min_radius, max_radius)

def generate_circles():
    centers = generate_circle_centers(n_circles)
    circles = []
    for center in centers:
        radius = generate_circle_radius()
        circle = center.buffer(radius)
        if polygon.contains(circle):
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

'''
# Plot the polygon
xp,yp = polygon.exterior.xy
plt.plot(xp,yp)
    #plt.show()

# generate random points
points = generate_circle_centers(polygon, 20)

# Plot the list of points
xs = [point.x for point in points]
ys = [point.y for point in points]
plt.scatter(xs, ys,color="red")
plt.show()
'''

# define requirements
n_circles = 20
min_radius = 0.005
max_radius = 0.01

xp,yp = polygon.exterior.xy
plt.plot(xp,yp)

circles = generate_circles()
for circle in circles:
    xp,yp = circle.exterior.xy
    plt.plot(xp,yp)

plt.show()