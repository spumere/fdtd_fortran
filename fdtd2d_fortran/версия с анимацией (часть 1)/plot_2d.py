import numpy as np
import matplotlib.pyplot as plt
import struct
from scipy.signal import find_peaks

def extremum_coords(A, prom=0.1, ax="x"):
    L = 0.9
    N = len(A)
    x = np.linspace(-L/2, L/2, N, endpoint=True)
    center = N // 2
    if ax == "x":
        A_1d = A[center, :]
    if ax == "y":
        A_1d = A[:, center]
    peaks_idx, _ = find_peaks(A_1d, distance=15, prominence=prom)
    peak_positions = x[peaks_idx]
    return peak_positions[peak_positions >= 0]

# --- Чтение бинарного файла ---
filename = 'fields2d.bin'

with open(filename, 'rb') as f:
    ie, jb, ib, je, nmax = struct.unpack('5q', f.read(40))       # 5 × int64
    dx, dt = struct.unpack('2d', f.read(16))                      # 2 × float64
    t_res = np.frombuffer(f.read(8 * nmax), dtype=np.float64)     # nt × float64

    ex_res = np.zeros((ie, jb, nmax), dtype=np.float64)
    ey_res = np.zeros((ib, je, nmax), dtype=np.float64)
    hz_res = np.zeros((ie, je, nmax), dtype=np.float64)

    for n in range(nmax):
        ex_res[:, :, n] = np.frombuffer(f.read(8 * ie * jb), dtype=np.float64).reshape((ie, jb), order='F')
        ey_res[:, :, n] = np.frombuffer(f.read(8 * ib * je), dtype=np.float64).reshape((ib, je), order='F')
        hz_res[:, :, n] = np.frombuffer(f.read(8 * ie * je), dtype=np.float64).reshape((ie, je), order='F')

# --- Два среза на одном плоте (2 строки x 3 столбца) ---
steps = [299, nmax - 1]
row_labels = [r'$\tau = \tau_1$', r'$\tau = \tau_2$']

fig, axes = plt.subplots(2, 3, figsize=(12, 10))

for row, (step, label) in enumerate(zip(steps, row_labels)):
    im_ex = axes[row, 0].imshow(ex_res[:, :, step].T, cmap='jet', vmin=-80, vmax=80,
                                origin='lower', interpolation='nearest', aspect='equal')
    im_ey = axes[row, 1].imshow(ey_res[:, :, step].T, cmap='jet', vmin=-80, vmax=80,
                                origin='lower', interpolation='nearest', aspect='equal')
    im_hz = axes[row, 2].imshow(hz_res[:, :, step].T, cmap='jet', vmin=-0.2, vmax=0.2,
                                origin='lower', interpolation='nearest', aspect='equal')

    plt.colorbar(im_ex, ax=axes[row, 0], shrink=0.6)
    plt.colorbar(im_ey, ax=axes[row, 1], shrink=0.6)
    plt.colorbar(im_hz, ax=axes[row, 2], shrink=0.6)

    for col in range(3):
        axes[row, col].set_xticks([])
        axes[row, col].set_yticks([])
        axes[row, col].set_aspect('equal')

    axes[row, 0].set_title(f'$E_x$, {label}')
    axes[row, 1].set_title(f'$E_y$, {label}')
    axes[row, 2].set_title(f'$H_z$, {label}')

plt.tight_layout()
plt.savefig('fields_combined.png', dpi=150)
plt.show()

# --- Расчёт расстояний между экстремумами для обоих шагов ---
for step, label in zip(steps, row_labels):
    ex_step = ex_res[:, :, step]
    ey_step = ey_res[:, :, step]
    hz_step = hz_res[:, :, step]

    dx_ex = np.mean(np.diff(extremum_coords(ex_step)))
    dx_ey = np.mean(np.diff(extremum_coords(ey_step, ax='y')))
    dx_hz = np.mean(np.diff(extremum_coords(hz_step)))

    print(f"--- {label} ---")
    print(f"  Среднее расстояние между максимумами Ex: {dx_ex:.8e} м")
    print(f"  Среднее расстояние между максимумами Ey: {dx_ey:.8e} м")
    print(f"  Среднее расстояние между максимумами Hz: {dx_hz:.8e} м")
