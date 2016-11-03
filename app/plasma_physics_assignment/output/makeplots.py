from math import pi, sin, cos, sqrt
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


DEFAULT_FIGSIZE = 6, 4
plt.rcParams.update({'font.size': 12,
                     'font.family': 'serif',
                     'text.usetex': True,
                     'text.latex.preamble': [r'\usepackage{gensymb}',
                                             r'\usepackage{lmodern}']})

# Solenoid B field
data = np.loadtxt('solenoid_B_onaxis.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('$z$ (m)')
plt.ylabel('$B_z$ (T)')
plt.xlim([-0.16, 0.17])
plt.tight_layout()
plt.savefig('solenoid_B_onaxis.pdf')

# Hexapole B field
data = np.loadtxt('hexapole_B.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.subplot(121)
plt.quiver(data[:, 0], data[:, 1], data[:, 2], data[:, 3])
plt.gca().set_aspect('equal')
plt.gca().locator_params(nbins=7)
plt.xlim([-0.039, 0.039])
plt.ylim([-0.039, 0.039])
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.subplot(122)
data.shape = 16, 16, 4
Bmid = midpoint_values(data[:, :, 2:4])
Bnorm = np.sqrt(Bmid[:, :, 0]**2 + Bmid[:, :, 1]**2)
plt.pcolormesh(data[:, :, 1], data[:, :, 0], Bnorm)
plt.gca().locator_params(nbins=7)
plt.xlim([-0.039, 0.039])
plt.ylim([-0.039, 0.039])
plt.colorbar(label=r'$|\mathbf{B}|$ (T)')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('hexapole_B.pdf')

# Total B field at (x, 0, z)
data = np.loadtxt('total_B_xz.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.subplot(211)
plt.quiver(data[:, 1], data[:, 0], data[:, 3], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlim([-0.16, 0.17])
plt.ylim([-0.039, 0.039])
plt.xlabel('$z$ (m)')
plt.ylabel('$x$ (m)')
plt.subplot(212)
data.shape = 21, 9, 4
Bmid = midpoint_values(data[:, :, 2:4])
Bnorm = np.sqrt(Bmid[:, :, 0]**2 + Bmid[:, :, 1]**2)
plt.pcolormesh(data[:, :, 1], data[:, :, 0], Bnorm)
plt.xlim([-0.16, 0.17])
plt.gca().set_aspect('equal')
plt.colorbar(label=r'$|\mathbf{B}|$ (T)')
plt.xlabel('z (m)')
plt.ylabel('x (m)')
plt.tight_layout()
plt.savefig('total_B_xz.pdf')

# Total B field at (0, y, z)
data = np.loadtxt('total_B_yz.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.subplot(211)
plt.quiver(data[:, 1], data[:, 0], data[:, 3], data[:, 2])
plt.gca().set_aspect('equal')
plt.xlim([-0.16, 0.17])
plt.ylim([-0.039, 0.039])
plt.xlabel('$z$ (m)')
plt.ylabel('$y$ (m)')
plt.subplot(212)
data.shape = 21, 9, 4
Bmid = midpoint_values(data[:, :, 2:4])
Bnorm = np.sqrt(Bmid[:, :, 0]**2 + Bmid[:, :, 1]**2)
plt.pcolormesh(data[:, :, 1], data[:, :, 0], Bnorm)
plt.xlim([-0.16, 0.17])
plt.colorbar(label=r'$|\mathbf{B}|$ (T)')
plt.gca().set_aspect('equal')
plt.xlabel('$z$ (m)')
plt.ylabel('$y$ (m)')
plt.tight_layout()
plt.savefig('total_B_yz.pdf')

# Single particle trajectory
data = np.loadtxt('particle_trajectory_10.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.plot(data[:, 1], data[:, 2])
#plt.gca().set_aspect('equal')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.savefig('particle_trajectory_10.pdf')
data = np.loadtxt('particle_trajectory_11.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.plot(data[:, 1], data[:, 2])
#plt.gca().set_aspect('equal')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.savefig('particle_trajectory_11.pdf')
data = np.loadtxt('particle_trajectory_12.csv', delimiter=';')
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.plot(data[:, 1], data[:, 2])
#plt.gca().set_aspect('equal')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.tight_layout()
plt.savefig('particle_trajectory_12.pdf')

data = np.loadtxt('z1_collision_points.csv', delimiter=';')
n, _ = data.shape
print('z1: ' + str(n))
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.hist2d(data[:, 0], data[:, 1], bins=50)
plt.colorbar(label='Number of collisions')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.tight_layout()
plt.savefig('z1_collision_points.pdf')

data = np.loadtxt('z2_collision_points.csv', delimiter=';')
n, _ = data.shape
print('z2: ' + str(n))
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.hist2d(data[:, 0], data[:, 1], bins=50)
plt.colorbar(label='Number of collisions')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.tight_layout()
plt.savefig('z2_collision_points.pdf')

data = np.loadtxt('cylinder_collision_points.csv', delimiter=';')
n, _ = data.shape
print('radial: ' + str(n))
plt.figure(figsize=DEFAULT_FIGSIZE)
plt.hist2d(np.arctan2(data[:, 1], data[:, 0]), data[:, 2], bins=100)
plt.colorbar(label='Number of collisions')
plt.xlabel(r'$\phi$ (rad)')
plt.ylabel('$z$ (m)')
plt.tight_layout()
plt.savefig('cylinder_collision_points.pdf')

# velocity component plot
nonlost = np.loadtxt('nonlost_velocity.csv', delimiter=';')[::30, :]
lost = np.loadtxt('lost_velocity.csv', delimiter=';')[::30, :]
plt.figure(figsize=(6,6))
pnonlost, = plt.plot(np.abs(nonlost[:, 0]), nonlost[:, 1], 'bs', label='Nonlost electrons', markersize=5, markeredgecolor='none')
plost, = plt.plot(np.abs(lost[:, 0]), lost[:, 1], 'ro', label='Lost electrons', markersize=3, markeredgecolor='none')

lower_angle = np.min(np.arctan(nonlost[:, 1] / np.abs(nonlost[:, 0])))
upper_angle = np.max(np.arctan(lost[:, 1] / np.abs(lost[:, 0])))
magnitude = sqrt(nonlost[0, 1]**2 + nonlost[0, 0]**2)

lower_line, = plt.plot([0, magnitude*cos(lower_angle)], [0, magnitude*sin(lower_angle)], label="{0:.1f}\\degree~line".format(lower_angle / pi * 180))
upper_line, = plt.plot([0, magnitude*cos(upper_angle)], [0, magnitude*sin(upper_angle)], label="{0:.1f}\\degree~line".format(upper_angle / pi * 180))

plt.xlabel(r'Radial velocity $v_\parallel$ (m/s)')
plt.ylabel(r'Transverse velocity $v_\perp$ (m/s)')
plt.legend(handles=[pnonlost, plost, lower_line, upper_line])
plt.tight_layout()
plt.savefig('loss_cone.pdf')

#plt.show()
