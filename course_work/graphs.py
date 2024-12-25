import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

with open('buf', 'r', encoding='utf-8') as buf_file:
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
    
    x = list(map(float, x))
    speed = list(map(float, speed))
    pressure = list(map(float, pressure))
    density = list(map(float, density))
    mach_number = list(map(float, mach_number))



fig, ax = plt.subplots(2, 2, figsize=(12, 12))

ax[0][0].plot(x, speed, label='Скорость')
ax[0][0].set_title('Скорость воздуха через сопло')
ax[0][0].set_xlabel('Продольная координата x, м', fontsize=14)
ax[0][0].set_ylabel('Cкорость u, м/c', fontsize=14)
# ax[0][0].set_yticks(np.append(ax[0][0].get_yticks(), [min(speed), max(speed)]))
ax[0][0].set_yticks(update_scale(ax[0][0].get_yticks(), min(speed), max(speed)))

ax[0][1].plot(x, mach_number, label='Число Маха')
ax[0][1].set_title('Число Маха')
ax[0][1].set_xlabel('Продольная координата x, м', fontsize=14)
ax[0][1].set_ylabel('Число Маха M', fontsize=14)
# ax[0][1].set_yticks(np.append(ax[0][1].get_yticks()[:-1], [min(mach_number), max(mach_number)]))
ax[0][1].set_yticks(update_scale(ax[0][1].get_yticks(), min(mach_number), max(mach_number)))

ax[1][0].plot(x, pressure, label='Давление')
ax[1][0].set_title('Давление воздуха')
ax[1][0].set_xlabel('Продольная координата x, м', fontsize=14)
ax[1][0].set_ylabel('Давление P, Па', fontsize=14)
# ax[1][0].set_yticks(np.append(ax[1][0].get_yticks()[:-1], [min(pressure), max(pressure)]))
ax[1][0].set_yticks(update_scale(ax[1][0].get_yticks(), min(pressure), max(pressure)))

ax[1][1].plot(x, density, label='Плотность')
ax[1][1].set_title('Плотность воздуха')
ax[1][1].set_xlabel('Продольная координата x, м', fontsize=14)
ax[1][1].set_ylabel('Плотность ρ, кг/м3', fontsize=14)
# ax[1][1].set_yticks(np.append(ax[1][1].get_yticks()[:-1], [min(density), max(density)]))
ax[1][1].set_yticks(update_scale(ax[1][1].get_yticks(), min(density), max(density)))

plt.subplots_adjust(hspace=0.3, wspace=0.3, top=0.97, bottom=0.05, left=0.08, right=0.98)
# plt.tight_layout()
plt.show()
