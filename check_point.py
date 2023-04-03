from shapely.geometry import Point, Polygon

# Define the polygon and the points to check
polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
points = [Point(0.5, 0.5), Point(2, 2)]

# Check if each point is inside the polygon
for point in points:
    if polygon.contains(point):
        print(f"{point} is inside the polygon")
    else:
        print(f"{point} is outside the polygon")


