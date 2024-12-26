import matplotlib.pyplot as plt
import numpy as np

def update_scale(scale, mn, mx):
    scale = np.delete(scale, np.argmax(scale))
    scale = np.delete(scale, np.argmin(np.abs(scale - mn)))
    scale = np.delete(scale, np.argmin(np.abs(scale - mx)))
    scale = np.append(scale, mn)
    scale = np.append(scale, mx)
    return scale

x = 1
speed = []
pressure = []
density = []
mach_number = []
r0 = 0
a = 0
x0 = 0
x1 = 0
radius = []

with open('buf', 'r', encoding='utf-8') as buf_file:
    tmp = buf_file.readline().split(' ')
    n = int(buf_file.readline())
    x = buf_file.readline().split(' ')
    speed = buf_file.readline().split(' ')
    pressure = buf_file.readline().split(' ')
    density = buf_file.readline().split(' ')
    mach_number = buf_file.readline().split(' ')
    
    x.pop()
    speed.pop()
    pressure.pop()
    density.pop()
    mach_number.pop()
    
    r0 = float(tmp[0])
    a = float(tmp[1])
    x = list(map(float, x))
    speed = list(map(float, speed))
    pressure = list(map(float, pressure))
    density = list(map(float, density))
    mach_number = list(map(float, mach_number))
    x0 = x[0]
    x1 = x[n-1]
    radius = [r0 + a*(t - x0) for t in x]

fig, ax = plt.subplots(2, 3, figsize=(18, 12))

ax[0][0].plot(x, radius, label='Профиль канала')
ax[0][0].set_title('Профиль канала')
ax[0][0].set_xlabel('Продольная координата x, м', fontsize=14)
ax[0][0].set_ylabel('Радиус r, м', fontsize=14)
ax[0][0].set_yticks(update_scale(ax[0][0].get_yticks(), min(radius), max(radius)))
ax[0][0].grid(color='black', linewidth=0.5)

ax[0][1].plot(x, speed, label='Скорость')
ax[0][1].set_title('Скорость воздуха')
ax[0][1].set_xlabel('Продольная координата x, м', fontsize=14)
ax[0][1].set_ylabel('Cкорость u, м/c', fontsize=14)
ax[0][1].set_yticks(update_scale(ax[0][1].get_yticks(), min(speed), max(speed)))
ax[0][1].grid(color='black', linewidth=0.5)

ax[0][2].plot(x, mach_number, label='Число Маха')
ax[0][2].set_title('Число Маха')
ax[0][2].set_xlabel('Продольная координата x, м', fontsize=14)
ax[0][2].set_ylabel('Число Маха M', fontsize=14)
ax[0][2].set_yticks(update_scale(ax[0][2].get_yticks(), min(mach_number), max(mach_number)))
ax[0][2].grid(color='black', linewidth=0.5)

ax[1][0].set_visible(False)

ax[1][1].plot(x, pressure, label='Давление')
ax[1][1].set_title('Давление воздуха')
ax[1][1].set_xlabel('Продольная координата x, м', fontsize=14)
ax[1][1].set_ylabel('Давление P, Па', fontsize=14)
ax[1][1].set_yticks(update_scale(ax[1][1].get_yticks(), min(pressure), max(pressure)))
ax[1][1].grid(color='black', linewidth=0.5)

ax[1][2].plot(x, density, label='Плотность')
ax[1][2].set_title('Плотность воздуха')
ax[1][2].set_xlabel('Продольная координата x, м', fontsize=14)
ax[1][2].set_ylabel('Плотность ρ, кг/м3', fontsize=14)
ax[1][2].set_yticks(update_scale(ax[1][2].get_yticks(), min(density), max(density)))
ax[1][2].grid(color='black', linewidth=0.5)

plt.subplots_adjust(hspace=0.3, wspace=0.3, top=0.97, bottom=0.05, left=0.08, right=0.98)
plt.show()
