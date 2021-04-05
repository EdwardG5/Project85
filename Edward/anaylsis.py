from collections import Counter


data = ""
with open("HCOV19-ENGLAND-081220-A.fasta", "r") as f:
    data = f.read()

data = data.split(">")
data = data[1:]

for i in range(len(data)):
    data[i] = data[i].split("8")
    data[i] = data[i][-1]
    data[i] = "".join(data[i].split())

def consensus(data):
    consensus = ""
    for i in range(len(data[0])):
        characters = list(map(lambda x: x[i], data))
        c = Counter(characters)
        mode = c.most_common(1)[0][0]
        consensus += mode
    return consensus

cS = consensus(data)

info = [() for _ in range(len(cS))]

for i in range(len(info)):
    info[i] += (cS[i],)
    disagrees = []
    for x in range(len(data)):
        if data[x][i] != cS[i]:
            disagrees.append(x)
    info[i] += (disagrees,)

numberWhichDisagree = list(map(lambda x: len(x[1]), info))

c = Counter(numberWhichDisagree)



####################################################################################

import numpy as np
import matplotlib.pyplot as plt


labels, values = zip(*sorted(c.items(), key=lambda x: x[0]))

fig = plt.figure()


# Make histogram of values with frequencies
plot1 = fig.add_subplot(121)
width = 1
plot1.bar(labels, values, width)

# Make a cdf
plot2 = fig.add_subplot(122)
l2 = (-1,)+labels # Add an extra value, otherwise it doesn't start at 0
v2 = (0,)+values
v2 = np.array(v2)
v2 = v2/len(cS) # Normalise
cy = np.cumsum(v2)
plot2.plot(l2, cy)

plt.show()