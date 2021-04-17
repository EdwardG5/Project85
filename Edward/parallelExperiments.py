from multiprocessing import Pool
import math

l = list(range(5))

def square(x):
    l = list(range(5000000))
    for i in range(len(l)):
        l[i] = i**2
    return sum(l)

# list(map(square, l))

# with Pool() as P:
#     P.map(square, l)

def f(x):
    print("f called")
    for i in range(10):
        x = x * x
    return x

import timeit

PARALLEL = 0

if __name__ == "__main__":
    print("Starting:\n")
    timeit.timeit(lambda: list(map(square, [1])), number=1)
    timeit.timeit(lambda: Pool().map(square, [1]), number=1)
    l = list(range(40))
    print("Sequential:")
    print(timeit.timeit(lambda: list(map(square, l)), number=1))
    print("\nParallel:")
    print(timeit.timeit(lambda: Pool().map(square, l), number=1))


"""
20:

Sequential:
25.736256281030364

Parallel:
17.674639272037894
"""

"""

"""