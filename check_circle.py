from shapely.geometry import Point, Polygon

# Define the polygon and the circle
polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
circle_center = Point(0.5, 0.5)
circle_radius = 0.25

# Check if the circle is entirely contained inside the polygon
circle = circle_center.buffer(circle_radius)
if polygon.contains(circle):
    print("The circle is entirely contained inside the polygon")
else:
    print("The circle is not entirely contained inside the polygon")