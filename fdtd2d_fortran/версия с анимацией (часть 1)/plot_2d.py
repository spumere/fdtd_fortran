import numpy as np
import matplotlib.pyplot as plt
import struct
from scipy.signal import find_peaks

def extremum_coords(A, prom = 0.1, ax = "x"):
    L = 0.9
    N = len(A)
    x = np.linspace(-L/2, L/2, N, endpoint = True)
    center = N // 2
    if ax == "x":
        A_1d = A[center, :]
    if ax == "y":
        A_1d = A[:, center]
    peaks_idx, _ = find_peaks(A_1d,  distance=15, prominence=prom)
    peak_positions = x[peaks_idx]
    return peak_positions[peak_positions >= 0]

# --- Чтение бинарного файла ---
filename = 'fields2d.bin'

with open(filename, 'rb') as f:
    # Заголовок: 5 целых по 4 байта
    ie, jb, ib, je, nmax = struct.unpack('5i', f.read(20))
    # dx, dt — 2 вещественных по 4 байта
    dx, dt = struct.unpack('2f', f.read(8))
    # Временной массив
    t_res = np.frombuffer(f.read(4 * nmax), dtype=np.float32)

    # Поля: для каждого шага ex(ie,jb), ey(ib,je), hz(ie,je)
    ex_res = np.zeros((ie, jb, nmax), dtype=np.float32)
    ey_res = np.zeros((ib, je, nmax), dtype=np.float32)
    hz_res = np.zeros((ie, je, nmax), dtype=np.float32)

    for n in range(nmax):
        ex_res[:, :, n] = np.frombuffer(f.read(4 * ie * jb), dtype=np.float32).reshape((ie, jb), order='F')
        ey_res[:, :, n] = np.frombuffer(f.read(4 * ib * je), dtype=np.float32).reshape((ib, je), order='F')
        hz_res[:, :, n] = np.frombuffer(f.read(4 * ie * je), dtype=np.float32).reshape((ie, je), order='F')

# --- Визуализация, как в MATLAB ---
fig, axes = plt.subplots(1, 3, figsize=(10, 6))

# Начальные изображения
im_ex = axes[0].imshow(ex_res[:, :, 0].T, cmap='jet', vmin=-80, vmax=80,
                       origin='lower', interpolation='nearest', aspect='equal')
im_ey = axes[1].imshow(ey_res[:, :, 0].T, cmap='jet', vmin=-80, vmax=80,
                       origin='lower', interpolation='nearest', aspect='equal')
im_hz = axes[2].imshow(hz_res[:, :, 0].T, cmap='jet', vmin=-0.2, vmax=0.2,
                       origin='lower', interpolation='nearest', aspect='equal')

# Цветовые шкалы
plt.colorbar(im_ex, ax=axes[0], shrink=0.4)
plt.colorbar(im_ey, ax=axes[1], shrink=0.4)
plt.colorbar(im_hz, ax=axes[2], shrink=0.4)

# Оформление осей
for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

axes[0].set_title('Ex at time step = 0')
axes[1].set_title('Ey at time step = 0')
axes[2].set_title('Hz at time step = 0')

plt.tight_layout()
plt.pause(0.001)

# --- Анимация ---
for n in range(nmax):
    im_ex.set_data(ex_res[:, :, n].T)
    im_ey.set_data(ey_res[:, :, n].T)
    im_hz.set_data(hz_res[:, :, n].T)

    axes[0].set_title(f'Ex at time step = {n+1}')
    axes[1].set_title(f'Ey at time step = {n+1}')
    axes[2].set_title(f'Hz at time step = {n+1}')

    plt.pause(0.001)

plt.ioff()
plt.show()
ex_res = ex_res[:, :, nmax-1]
ey_res = ey_res[:, :, nmax-1]
hz_res = hz_res[:, :, nmax-1]

print(f"Среднее расстояние между экстремумами функции Ex: {np.mean(np.diff(extremum_coords(ex_res))):.8e} м")
print(f"Среднее расстояние между экстремумами функции Ey: {np.mean(np.diff(extremum_coords(ey_res, ax = "y"))):.8e} м")
print(f"Среднее расстояние между экстремумами функции Hz: {np.mean(np.diff(extremum_coords(hz_res))):.8e} м")