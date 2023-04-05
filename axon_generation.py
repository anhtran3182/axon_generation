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

def get_axon_info(axons):
    center_x = []
    center_y = []
    radius = []
    for axon in axons:
        center_x.append(axon.centroid.x)
        center_y.append(axon.centroid.y)
        radius.append(axon.exterior.coords[0][0] - axon.centroid.x)
    return center_x, center_y, radius
    
def save_axon_info(axons, filename):
    center_x, center_y, radius = get_axon_info(axons)

    with open(filename, mode='w', newline='') as csvfile:
        fieldnames = ['contour_num', 'centroid_x', 'centroid_y', 'centroid_z', 'circular_diameter']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(len(axons)):
            writer.writerow({'contour_num': i+1, 'centroid_x': center_x[i], 'centroid_y': center_y[i], 'centroid_z': 0, 'circular_diameter': 2*radius[i]})

def plot():
    # plot the shape and axons
    xp,yp = shape.exterior.xy
    plt.plot(xp,yp)

    for each_axon in axons:
        xp,yp = each_axon.exterior.xy
        plt.plot(xp,yp)

    plt.show()


# make the polygon 
boundary = generate_boundary('Results_Sample22209-20x.csv')
shape = Polygon(boundary)

# define requirements
axons_num = 100
min_radius = 1
max_radius = 5
axons = []

# redefine the shape so that the generated circles do not touch the boundary 
buffer_distance = min_radius * 0.5
shape = shape.buffer(-buffer_distance)

axons = generate_circles()
save_axon_info(axons, 'axon_info.csv')

plot()

'''
Axon outputs:
   contour_num
   centroid_x
   centroid_y
   centroid_z = all zeros
   circular diameter


need a box of saline
'''
