import matplotlib.pyplot as plt

accuracy = [x for x in range(0, 16)]
rotation_iterations = [2, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8]
power_iterations = [7, 14, 21, 28, 34, 41, 48, 55, 62, 69, 76, 83, 90, 97, 103, 110]

plt.plot(accuracy, rotation_iterations, label="Метод вращений")
plt.plot(accuracy, power_iterations, label="Степенной метод")

plt.xlabel('Точность, e^(-x)', fontsize=14)
plt.ylabel('Количество итераций', fontsize=14)

plt.legend()

plt.show()