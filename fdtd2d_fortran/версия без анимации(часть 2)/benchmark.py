# Бенчмарк векторизации FDTD 2D.

# Что делает:
#   1. Компилирует две версии: scalar (без векторизации) и vector (с векторизацией).
#   2. Прогоняет обе для списка nmax.
#   3. Считает S = t_scalar / t_vector.
#   4. Пишет results/benchmark.csv.
#   5. Строит results/benchmark.png с двумя графиками.
#   6. Для последнего прогона сохраняет:
#        - results/metadata.txt   (λ, t, расстояния между экстремумами)
#        - results/fields2d_final.png (картинка последнего кадра)

import os
import csv
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import struct

SOURCE = 'fdtd2d.f90'
RESULTS_DIR = 'results'
EXE_SCALAR = 'fdtd2d_scalar.exe'
EXE_VECTOR = 'fdtd2d_vector.exe'

NMAX_LIST = np.arange(300, 5100, 200)

REPEATS = 3

C = 2.99792458e8
dx = 3.0e-3
dt = dx/(2.0*C) * 1e9

# Флаги компиляции
# Cравниваем только эффект векторизации, всё остальное одинаково.
# -O2 без векторизации vs -O2 с векторизацией.
FLAGS_COMMON = ['-O2', '-march=native']
FLAGS_SCALAR = FLAGS_COMMON + ['-fno-tree-vectorize', '-fno-tree-slp-vectorize']
FLAGS_VECTOR = FLAGS_COMMON + ['-ftree-vectorize', '-funroll-loops']

def compile_exe(flags, source, out_exe):
    # Компилирует source с заданными флагами в out_exe.
    cmd = ['gfortran'] + flags + [source, '-o', out_exe]
    print(f"  {' '.join(cmd)}")
    subprocess.run(cmd, capture_output=True, text=True)

def run_once(exe, nmax):
    # Запускает exe nmax, возвращает время в секундах.
    r = subprocess.run([exe, str(nmax)], capture_output=True, text=True)
    # Ищем строку "Calculation time = 0.1234 s"
    for line in r.stdout.splitlines():
        if 'Calculation time' in line:
            # "Calculation time =  0.1234  s"
            t_str = line.split('=')[1].split('s')[0].strip()
            return float(t_str)

def run_best_of(exe, nmax, repeats=REPEATS):
    # Запускает REPEATS раз, возвращает минимальное время
    times = [run_once(exe, nmax) for _ in range(repeats)]
    return min(times)

def extremum_coords(A, prom=0.1, ax="x"):
    L = 0.9
    N = A.shape[0]
    x = np.linspace(-L / 2, L / 2, N, endpoint=True)
    center = N // 2
    if ax == "x":
        A_1d = A[center, :]
    else:
        A_1d = A[:, center]
    peaks_idx, _ = find_peaks(A_1d, distance=15, prominence=prom)
    peak_positions = x[peaks_idx]
    return peak_positions[peak_positions >= 0]


def read_bin(filename):
    # Читает fields2d.bin с одним кадром.
    with open(filename, 'rb') as f:
        ie, jb, ib, je, nmax = struct.unpack('5i', f.read(20))
        dx, dt = struct.unpack('2f', f.read(8))
        t_final = struct.unpack('f', f.read(4))[0]
        ex = np.frombuffer(f.read(4 * ie * jb), dtype=np.float32).reshape((ie, jb), order='F')
        ey = np.frombuffer(f.read(4 * ib * je), dtype=np.float32).reshape((ib, je), order='F')
        hz = np.frombuffer(f.read(4 * ie * je), dtype=np.float32).reshape((ie, je), order='F')
    return dict(ie=ie, jb=jb, ib=ib, je=je, dx=dx, dt=dt,
                t_final=t_final, ex=ex, ey=ey, hz=hz)


def save_metadata(data, nmax, path):
    # Считает экстремумы и пишет metadata.txt.
    ex, ey, hz = data['ex'], data['ey'], data['hz']

    d_ex = np.mean(np.diff(extremum_coords(ex)))
    d_ey = np.mean(np.diff(extremum_coords(ey, ax='y')))
    d_hz = np.mean(np.diff(extremum_coords(hz)))

    lambda_m = C / 5.0e9

    with open(path, 'w', encoding='utf-8') as f:
        f.write(f"nmax = {nmax}\n")
        f.write(f"ie = {data['ie']}\n")
        f.write(f"je = {data['je']}\n")
        f.write(f"dx = {data['dx']:.6e} м\n")
        f.write(f"dt = {data['dt']:.6e} с\n")
        f.write(f"lambda = {lambda_m:.8e} м\n")
        f.write(f"t = {data['t_final']:.6f} нс\n")
        f.write(f"mean_distance_ex = {d_ex:.8e} м\n")
        f.write(f"mean_distance_ey = {d_ey:.8e} м\n")
        f.write(f"mean_distance_hz = {d_hz:.8e} м\n")
    print(f"Информация о последнем изображении для nmax -> {path}")


def save_final_image(data, path):
    # Сохраняет PNG с последним кадром.
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))

    im_ex = axes[0].imshow(data['ex'].T, cmap='jet', vmin=-80, vmax=80,
                           origin='lower', interpolation='nearest', aspect='equal')
    im_ey = axes[1].imshow(data['ey'].T, cmap='jet', vmin=-80, vmax=80,
                           origin='lower', interpolation='nearest', aspect='equal')
    im_hz = axes[2].imshow(data['hz'].T, cmap='jet', vmin=-0.2, vmax=0.2,
                           origin='lower', interpolation='nearest', aspect='equal')

    plt.colorbar(im_ex, ax=axes[0], shrink=0.7)
    plt.colorbar(im_ey, ax=axes[1], shrink=0.7)
    plt.colorbar(im_hz, ax=axes[2], shrink=0.7)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

    axes[0].set_title('Ex')
    axes[1].set_title('Ey')
    axes[2].set_title('Hz')

    fig.suptitle(f"Последнее распределение, t = {data['t_final']:.4f} нс", fontsize=14)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close(fig)

# ============================ ГРАФИК ============================

def plot_benchmark(rows, path):
    # Строит два графика: времена и S = t_scalar / t_vector.
    nmax = [r['nmax']*dt for r in rows]
    t_s = [r['t_scalar'] for r in rows]
    t_v = [r['t_vector'] for r in rows]
    S = [r['S'] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Левый: времена
    axes[0].plot(nmax, t_s, 'o-', label='scalar (-O2 -fno-tree-vectorize)', color='C3')
    axes[0].plot(nmax, t_v, 's-', label='vector (-O2 -ftree-vectorize)', color='C0')
    axes[0].set_xlabel('Моделируемый временной промежуток, нс')
    axes[0].set_ylabel('Время расчёта, с')
    axes[0].set_title('Время расчёта')
    axes[0].legend()
    axes[0].grid(True, alpha=0.4)

    # Правый: ускорение
    axes[1].plot(nmax, S, 'D-', color='C2')
    axes[1].axhline(1.0, color='gray', linestyle='--',
                    label='S = 1')
    axes[1].set_xlabel('Моделируемый временной промежуток, нс')
    axes[1].set_ylabel('S = $t_{vector} / t_{scalar}$')
    axes[1].set_title('Ускорение от векторизации')
    axes[1].legend()
    axes[1].grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close(fig)
    print(f"график -> {path}")



os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 70)
print("Компиляция")
print("=" * 70)
print("Скалярный вариант:")
compile_exe(FLAGS_SCALAR, SOURCE, EXE_SCALAR)
print("Векторный вариант:")
compile_exe(FLAGS_VECTOR, SOURCE, EXE_VECTOR)

print()
print("=" * 70)
print(f"Выполнение для t = {NMAX_LIST * dt} нс")
print(f"Повторов на каждый запуск: {REPEATS} (результат - минимальное время)")
print("=" * 70)

rows = []
for nmax in NMAX_LIST:
    print(f"\n---t = {(nmax * dt):.6f} нс---")
    print("\tСкалярный вариант: ")
    t_scalar = run_best_of(EXE_SCALAR, nmax)
    print(f"\tt_scalar = {t_scalar:.6f} с")
    print("\tВекторный вариант:")
    t_vector = run_best_of(EXE_VECTOR, nmax)
    print(f"\tt_vector = {t_vector:.6f} с")

    S = t_scalar / t_vector if t_vector > 0 else float('inf')
    print(f"\tS = {S:.4f}")

    rows.append(dict(nmax=nmax, t_scalar=t_scalar, t_vector=t_vector, S=S))

# CSV
csv_path = os.path.join(RESULTS_DIR, 'benchmark.csv')
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['nmax', 't_scalar', 't_vector', 'S'])
    writer.writeheader()
    writer.writerows(rows)
print(f"\nCSV -> {csv_path}")

# График
plot_benchmark(rows, os.path.join(RESULTS_DIR, 'benchmark.png'))

nmax_last = NMAX_LIST[-1]
run_once(EXE_VECTOR, nmax_last)

bin_path = os.path.join(RESULTS_DIR, 'fields2d.bin')

data = read_bin(bin_path)
save_metadata(data, nmax_last, os.path.join(RESULTS_DIR, 'metadata.txt'))
save_final_image(data, os.path.join(RESULTS_DIR, 'fields2d_final.png'))

