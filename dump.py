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


# Function for filtering data points
def moving_average_filter(N: int, n: int, f):
    # Declare an x array of length n, with steps of 1
    xn = np.arange(0, n)
    yn = []

    # For each step in the x array, calculate the moving average based on some value N
    for i in range(len(xn)):
        # Use a running sum to calculate the current data point
        curr = 0

        # Take the current point and add it to 7 previous points
        # If the index is negative then add 0 (i.e., don't do anything)
        for j in range(i, i - N + 1, -1):
            if j >= 0:
                curr += f(xn[j])

        # Complete the formula by multiplying the current point by 1/N and add it to the array
        curr *= 1/N
        yn.append(curr)

    return yn


# Since we have an 8 point moving average filter, we use a const N = 8
N = 8

# Declare an x array of length 60 (advice taken from hints), with steps of 1 for option A
xn_a = np.arange(0, 60)
yn_a_unsmoothed = [f(xn_a[i]) for i in range(len(xn_a))]
yn_a = moving_average_filter(N, len(xn_a), f)

# Plot data to figure 1
plt.figure(1)
plt.grid("on")
plt.plot(xn_a, yn_a)
plt.plot(xn_a, yn_a_unsmoothed)
plt.legend(["Filtered", "Original"])
plt.title("Question A")

# Derive data for option B and plot to figure 2
xn_b = np.arange(0, 90)
yn_b_unsmoothed = [g(xn_b[i]) for i in range(len(xn_b))]
yn_b = moving_average_filter(N, len(xn_b), g)

plt.figure(2)
plt.grid("on")
plt.title("Question B")
plt.plot(xn_b, yn_b)
plt.plot(xn_b, yn_b_unsmoothed)
plt.legend(["Filtered", "Original"])

# Show plots
plt.show()
