# Stochastic Events Course

import numpy as np
import math
import matplotlib.pyplot as plt


#HW3
# def randwalk(n):
#     x = 0
#     y = 0
#     time = [x]
#     position = [y]
#     for i in range(1, n + 1):
#         move = np.random.uniform(0, 1)
#         if move < 0.75:
#             x += 1
#             y += 1
#         if move > 0.75:
#             x += 1
#             y += -1
#         time.append(x)
#         position.append(y)
#
#     return [time, position]
# n = 10000
# p_array = np.zeros(n)
# for i in range(0, len(p_array)):
#     rw = randwalk(20)
#     p_array[i] = (rw[1][-1])



# plt.hist(p_array,bins=40,density=True)
# plt.xlabel("step")
# plt.ylabel("P")
# plt.title("p=0.75")
# plt.show()

#HW 4
# def f(t):
#     r=1
#     return (2/(math.pi*r**2))* np.sqrt(r**2 - t**2)
#
# t2 = np.arange(-5, 5.0, 0.02)
#
# plt.figure()
# plt.plot(t2, f(t2), 'k')
# plt.show()



# myMatrix = np.random.normal(size=(5000, 5000))
# for i in range(5000):
#     for j in range(5000):
#         myMatrix[j][i] = myMatrix[i][j]
#
# w , v = np.linalg.eig(myMatrix)
#
# plt.hist(w/(2*math.sqrt(5000)), bins=40, density=True)
# plt.plot(t2, f(t2), 'k')
# plt.legend("Analytical Result")
# plt.title("Histogram")
# plt.show()

# HW 5

n=20
l = np.sqrt(3)
sigma_s = 1
a = np.sqrt(2)
b = 1
########## DISTRIBUTION FUCNTIONS ############
# x = np.arange(-l,l,0.01)   # start,stop,step
# y = 1/2*l
# z = (1/(np.sqrt(2*np.pi*sigma_s))) * np.exp(-(x**2/2*sigma_s))
# w = (a/2) * np.exp(-a * abs(x))
# m = (b/ np.pi) * (1/(b**2 + x**2))
#
# plt.axhline(y=y, color='r', linestyle='-')
# plt.plot(x, z, x, w, x, m)
# plt.xlabel('x')  # string must be enclosed with quotes '  '
# plt.ylabel('y')
# plt.title('Distribution functions')
# plt.legend(['Normal', 'Gaussian', 'Laplace', 'Cauchy'])
# plt.show()
########## CHARACRETRISTIC FUCNTIONS ############
k = np.arange(-l*2, l*2, 0.01)   # start,stop,step
y = (np.sin(l*k) / (l*k))**n
z = (np.exp(-sigma_s*k**2/2))**n
w = ((a**2 / (a**2+k**2)))**n
m = (np.exp(-b * abs(k)))**n

# plt.axhline(y=y, color='r', linestyle='-')
plt.plot(k, y, k, z, k, w, k, m)
plt.xlabel('k')  # string must be enclosed with quotes '  '
plt.ylabel('P_3')
plt.title('Characteritic functions')
plt.legend(['sin(l*k) / (l*k)', 'exp(-sigma_s*k**2/2)', '(a**2 / (a**2+k**2))', 'exp(-b * abs(k))'])
plt.show()












