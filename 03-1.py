import numpy as np
import matplotlib.pyplot as plt

def u(x, e):
    if 0 <= x <= 100:
        return 100 * x / e
    elif 100 <= x <= 180:
        return 10000 / e + 4000 / e - 4000 / (e * (1 + (x - 100) / 40))

def sigma(x):
    if 0 <= x <= 100:
        return 100
    elif 100 <= x <= 180:
        return 100 / (1 + (x - 100) / 40)

if __name__ == '__main__':
    x = np.linspace(0, 180, 200)
    y = [sigma(xi) for xi in x]
    plt.plot(x, y)
    plt.ylabel('$\sigma$(cm)')
    plt.xlabel('$x$(cm)')
    plt.show()
