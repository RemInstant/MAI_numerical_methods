import matplotlib.pyplot as plt

accuracy = [x for x in range(0, 16)]
iterations = [5, 6, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, 26, 28, 29]
seidel = [2, 3, 5, 5, 6, 8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18]

plt.plot(accuracy, iterations, label="Метод итераций")
plt.plot(accuracy, seidel, label="Метод Зейделя")

plt.xlabel('Точность, e^(-x)', fontsize=14)
plt.ylabel('Количество итераций', fontsize=14)

plt.legend()

plt.show()