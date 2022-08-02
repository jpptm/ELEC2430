import numpy as np
import matplotlib.pyplot as plt

# Exercise 1 - 8 point Moving average filter

# Option A: x[n] = 2sin(pi*n/10 - pi/3) 0 <= n <= 40
# Option B: x[n] = 2cos(pi*n/10 - pi/3) 0 <= n <= 60

# Output generalisation
# y[n] = 1/N * sigma(x[n - i], i=[0, N-1])

# Option A's function


def f(n):
    return 2*np.sin(np.pi*n/10 - np.pi/3)

# Option B's function


def g(n):
    return 2*np.cos(np.pi*n/10 - np.pi/3)


N = 8
xn = np.arange(0, 90)
yn_unsmoothed = [f(xn[i]) for i in range(len(xn))]
yn = []

for i in range(len(xn)):

    curr = 0

    for j in range(i, i - N + 1, -1):
        if j >= 0:
            curr += f(xn[j])

    curr *= 1/N
    yn.append(curr)

plt.grid("on")
plt.plot(xn, yn)
plt.plot(xn, yn_unsmoothed)
plt.legend(["Filtered", "Original"])
plt.title("Question A")
plt.show()


def moving_average_filter(N: int, n: int, f):
    xn = np.arange(0, n)
    yn = []

    for i in range(len(xn)):

        curr = 0

        for j in range(i, i - N + 1, -1):
            if j >= 0:
                curr += f(xn[j])

        curr *= 1/N
        yn.append(curr)

    return yn


xn_b = np.arange(0, 90)
yn_b_unfiltered = [g(xn_b[i]) for i in range(len(xn_b))]
yn_b = moving_average_filter(N, len(xn_b), g)
plt.grid("on")
plt.legend(["Filtered", "Original"])
plt.title("Question B")
plt.plot(xn_b, yn_b)
plt.plot(xn_b, yn_b_unfiltered)
plt.show()
