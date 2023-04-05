import matplotlib.pyplot as plt
import numpy as np

def generate_center_position(center, radius):
    while True:
        x = np.random.uniform(-radius, radius)
        y = np.random.uniform(-radius, radius)
        
        if np.linalg.norm(np.array([x, y]) - center) < radius:
            # The point is within the arbitrary circle
            return np.array([x, y])

def is_circle_overlapping(center, radius, circles):
    for circle in circles:
        if np.linalg.norm(center - circle[:2]) < radius + circle[2]:
            return True
    return False

def generate_circles(center, radius, n_circles, min_radius, max_radius):
    circles = []

    # Generate the circles
    while len(circles) < n_circles:
        generated_radius = np.random.uniform(min_radius, max_radius)
        generated_center = generate_center_position(center, radius - generated_radius)

        # Check if the circle is within bounds and does not overlap with any existing circle
        if not is_circle_overlapping(generated_center, generated_radius, circles):
            circles.append(np.concatenate((generated_center, [generated_radius])))

    return circles

def main():
    # Generate the circles
    center = np.array([0, 0])
    radius = 5
    n_circles = 20
    min_radius = 0.5
    max_radius = 1
    circles = generate_circles(center, radius, n_circles, min_radius, max_radius)

    # Plot the circles
    fig, ax = plt.subplots()
    circle1 = plt.Circle((0, 0), radius, fill=False)
    ax.add_artist(circle1)
    for circle in circles:
        circle2 = plt.Circle(circle[:2], circle[2], color='blue', alpha=0.5)
        ax.add_artist(circle2)
        
    plt.xlim(-radius, radius)
    plt.ylim(-radius, radius)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

if __name__ == '__main__':
    main()
