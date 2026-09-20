import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

with open('vectorized_results/fields.bin', 'rb') as f:
    ih   = np.fromfile(f, dtype=np.int32, count=1)[0]
    nmax = np.fromfile(f, dtype=np.int32, count=1)[0]

    x     = np.fromfile(f, dtype=np.float32, count=ih)
    t     = np.fromfile(f, dtype=np.float32, count=nmax)
    ez    = np.fromfile(f, dtype=np.float32, count=nmax*ih).reshape(nmax, ih)
    ez_an = np.fromfile(f, dtype=np.float32, count=nmax*ih).reshape(nmax, ih)
    hy    = np.fromfile(f, dtype=np.float32, count=nmax*ih).reshape(nmax, ih)
    hy_an = np.fromfile(f, dtype=np.float32, count=nmax*ih).reshape(nmax, ih)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

def update(k):
    ax1.clear(); ax2.clear()

    ax1.plot(x, ez[k, :],    'r-',  label='Одномерный метод Йи')
    ax1.plot(x, ez_an[k, :], 'b--', label='Аналитическое решение')
    ax1.set_ylabel('$E_z$')
    ax1.set_ylim(-1, 1)
    ax1.set_title(f't = {t[k]:.2f} нс')
    ax1.legend(); ax1.grid(True)

    ax2.plot(x, hy[k, :],    'g-',  label='Одномерный метод Йи')
    ax2.plot(x, hy_an[k, :], 'm--', label= 'Аналитическое решение')
    ax2.set_ylabel('$H_y$')
    ax2.set_ylim(-3e-3, 3e-3)
    ax2.set_xlabel('x (м)')
    ax2.legend(); ax2.grid(True)

ani = FuncAnimation(fig, update, frames=int(nmax), interval=50, repeat=True)
plt.show()