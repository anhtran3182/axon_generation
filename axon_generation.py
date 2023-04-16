import csv
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, mapping

'''
This file: 'Results_Sample22209-20x.csv'

boundary_x is array for column 1
boundary_y is array for column 2
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
            if 'X' not in row[1] and 'Y' not in row[2]:
                boundary_x.append(row[1])  
                boundary_y.append(row[2])
    return list(zip(boundary_x, boundary_y))

def get_axon_info():
    axon_info = []
    pixels_per_micron = 2401.336 / 1000
    
    for i, axon in enumerate(axons):
        info = {
            "Contour_num": i+1,
            "Centroid_x": axon.centroid.x / pixels_per_micron,
            "Centroid_y": axon.centroid.y / pixels_per_micron,
            "Centroid_z": 0,
            "Circular Diameter (Area)": axon.exterior.coords[0][0] - axon.centroid.x
        }
        axon_info.append(info)
        
    return axon_info

def save_axon_info(filename):
    axon_info = get_axon_info()
    
    with open(filename, "w", newline="") as csvfile:
        fieldnames = ["Contour_num", "Centroid_x", "Centroid_y", "Centroid_z", "Circular Diameter (Area)"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for info in axon_info:
            writer.writerow(info)
            
def print_axon_info():
    axon_info = get_axon_info()
    print("Contour_num | Centroid X | Centroid Y | Diameter")
    for info in axon_info:
        print(f"{info['Contour_num']:9d} | {info['Centroid_x']:10.2f} | {info['Centroid_y']:10.2f} | {info['Circular Diameter (Area)']:8.2f}")

            
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

save_axon_info('axon_info.csv')
print_axon_info()
plot()

