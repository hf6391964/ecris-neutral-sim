import numpy as np
from matplotlib import pyplot as plt


def midpoint_values(arr):
    sy, sx, _ = arr.shape
    return 0.25 * (
        arr[0:-1, 0:-1, :] +
        arr[0:-1, 1:, :] +
        arr[1:, 0:-1, :] +
        arr[1:, 1:, :]
    )


# Solenoid B field
data = np.loadtxt('solenoid_B_onaxis.csv', delimiter=';')
plt.figure()
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('z (m)')
plt.ylabel('B (T)')
plt.xlim([-0.16, 0.17])
plt.savefig('solenoid_B_onaxis.pdf')

# Hexapole B field
data = np.loadtxt('hexapole_B.csv', delimiter=';')
plt.figure()
plt.subplot(121)
plt.quiver(data[:, 0], data[:, 1], data[:, 2], data[:, 3])
plt.gca().set_aspect('equal')
plt.gca().locator_params(nbins=7)
plt.xlim([-0.039, 0.039])
plt.ylim([-0.039, 0.039])
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.subplot(122)
data.shape = 16, 16, 4
Bmid = midpoint_values(data[:, :, 2:4])
Bnorm = np.sqrt(Bmid[:, :, 0]**2 + Bmid[:, :, 1]**2)
plt.pcolormesh(data[:, :, 1], data[:, :, 0], Bnorm)
plt.gca().locator_params(nbins=7)
plt.xlim([-0.039, 0.039])
plt.ylim([-0.039, 0.039])
plt.colorbar()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('hexapole_B.pdf')

# Total B field at (x, 0, z)
data = np.loadtxt('total_B_xz.csv', delimiter=';')
plt.figure()
plt.subplot(211)
plt.quiver(data[:, 1], data[:, 0], data[:, 3], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlim([-0.16, 0.17])
plt.ylim([-0.039, 0.039])
plt.xlabel('z (m)')
plt.ylabel('x (m)')
plt.subplot(212)
data.shape = 21, 9, 4
Bmid = midpoint_values(data[:, :, 2:4])
Bnorm = np.sqrt(Bmid[:, :, 0]**2 + Bmid[:, :, 1]**2)
plt.pcolormesh(data[:, :, 1], data[:, :, 0], Bnorm)
plt.xlim([-0.16, 0.17])
plt.gca().set_aspect('equal')
plt.colorbar()
plt.xlabel('z (m)')
plt.ylabel('x (m)')
plt.savefig('total_B_xz.pdf')

# Total B field at (0, y, z)
data = np.loadtxt('total_B_yz.csv', delimiter=';')
plt.figure()
plt.subplot(211)
plt.quiver(data[:, 1], data[:, 0], data[:, 3], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlim([-0.16, 0.17])
plt.ylim([-0.039, 0.039])
plt.xlabel('z (m)')
plt.ylabel('y (m)')
plt.subplot(212)
data.shape = 21, 9, 4
Bmid = midpoint_values(data[:, :, 2:4])
Bnorm = np.sqrt(Bmid[:, :, 0]**2 + Bmid[:, :, 1]**2)
plt.pcolormesh(data[:, :, 1], data[:, :, 0], Bnorm)
plt.xlim([-0.16, 0.17])
plt.colorbar()
plt.gca().set_aspect('equal')
plt.xlabel('z (m)')
plt.ylabel('y (m)')
plt.savefig('total_B_yz.pdf')

# Single particle trajectory
data = np.loadtxt('particle_trajectory_10.csv', delimiter=';')
plt.figure()
plt.subplot(121)
plt.plot(data[:, 1], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.subplot(122)
print(data[:, 0].shape)
plt.plot(data[:, 0] * 1e9, data[:, 4])
# plt.gca().set_aspect('equal')
plt.xlabel('t (ns)')
plt.ylabel('E (eV)')
plt.savefig('particle_trajectory_10.pdf')

data = np.loadtxt('z1_collision_points.csv', delimiter=';')
plt.figure()
plt.hist2d(data[:, 0], data[:, 1], bins=50)
plt.colorbar()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.savefig('z1_collision_points.pdf')

data = np.loadtxt('z2_collision_points.csv', delimiter=';')
plt.figure()
plt.hist2d(data[:, 0], data[:, 1], bins=50)
plt.colorbar()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.savefig('z2_collision_points.pdf')

data = np.loadtxt('cylinder_collision_points.csv', delimiter=';')
plt.figure()
plt.hist2d(np.arctan2(data[:, 1], data[:, 0]), data[:, 2], bins=100)
plt.colorbar()
plt.xlabel('phi (rad)')
plt.ylabel('z (m)')
plt.savefig('cylinder_collision_points.pdf')

plt.show()
