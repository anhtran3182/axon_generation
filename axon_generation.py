import csv
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, mapping

'''
This file: 'Results_Sample22209-20x.csv'

boundary_x is array for column 6
boundary_y is array for column 7
'''

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, mapping

def generate_circle_centers():
    minx, miny, maxx, maxy = shape.bounds
    while True:
        point = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
        if shape.contains(point):
            return point

def generate_circle_radius():
    return np.random.uniform(min_radius, max_radius)

def generate_circles():
    axons = []

    while len(axons) < axons_num:
        center = generate_circle_centers()
        radius = generate_circle_radius()
        new_axon = center.buffer(radius)
        if not any(existing_axon.intersects(new_axon) for existing_axon in axons) and shape.contains(new_axon):
            axons.append(new_axon)

    return axons

def generate_boundary(filename):
    boundary_x = []  
    boundary_y = []

    with open(filename, 'r') as file:
        reader = csv.reader(file)

        for row in reader:
            if 'X' not in row[5] and 'Y' not in row[6]:
                boundary_x.append(row[5])  
                boundary_y.append(row[6])
    return list(zip(boundary_x, boundary_y))
    
    
# make the polygon 
boundary = generate_boundary('Results_Sample22209-20x.csv')
shape = Polygon(boundary)

# define requirements
axons_num = 100
min_radius = 1
max_radius = 5

# redefine the shape so that the generated circles do not touch the boundary 
buffer_distance = min_radius * 0.5
shape = shape.buffer(-buffer_distance)

# plot the shape and axons
xp,yp = shape.exterior.xy
plt.plot(xp,yp)

circles = generate_circles()
for circle in circles:
    xp,yp = circle.exterior.xy
    plt.plot(xp,yp)

plt.show()