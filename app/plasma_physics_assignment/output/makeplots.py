import numpy as np
from matplotlib import pyplot as plt


# Solenoid B field
data = np.loadtxt('solenoid_B_onaxis.csv', delimiter=';')
plt.figure()
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('z (m)')
plt.ylabel('B (T)')
plt.xlim([-0.16, 0.17])
plt.savefig('solenoid_B_onaxis.png')

# Hexapole B field
data = np.loadtxt('hexapole_B.csv', delimiter=';')
plt.figure()
plt.quiver(data[:, 0], data[:, 1], data[:, 2], data[:, 3])
plt.gca().set_aspect('equal')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.savefig('hexapole_B.png')

# Total B field at (x, 0, z)
data = np.loadtxt('total_B_xz.csv', delimiter=';')
plt.figure()
plt.quiver(data[:, 1], data[:, 0], data[:, 3], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlim([-0.16, 0.17])
plt.xlabel('z (m)')
plt.ylabel('x (m)')
plt.savefig('total_B_xz.png')

# Total B field at (0, y, z)
data = np.loadtxt('total_B_yz.csv', delimiter=';')
plt.figure()
plt.quiver(data[:, 1], data[:, 0], data[:, 3], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlim([-0.16, 0.17])
plt.xlabel('z (m)')
plt.ylabel('y (m)')
plt.savefig('total_B_yz.png')

plt.show()
