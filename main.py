import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import numpy as np

dr = 0.01
#R = 1
#n_r = int(R / dr)
dt = 0.01
#total_time = 10
#n_t = int(total_time / dt)

num_steps = 100
l = np.empty(num_steps + 2, dtype=np.int16)

def fick(c, scheme):
    diff = 0.0001
    if scheme == "explicit":
        c[1:-1] = (diff * dt / dr**2) * (c[2:] - 2 * c[1:-1] + c[:-2]) + c[1:-1]
        c[0] = c[1]
        c[-1] = 0

    elif scheme == "implicit":
        alfa_i = np.zeros(num_steps + 2)
        betta_i = np.zeros(num_steps + 2)
        alfa_i[0] = 1
        betta_i[0] = 0

        a = -diff * dt / dr**2
        b = 1 + 2 * diff * dt / dr**2
        c_koef = -diff * dt / dr**2

        for i in range(1, len(l)):
            alfa_i[i] = (-a) / (b + c_koef * alfa_i[i - 1])
            betta_i[i] = (c[i] - c_koef * betta_i[i - 1]) / (b + c_koef * alfa_i[i - 1])
        alfa_i[-1] = 0
        betta_i[-1] = 0
        c[1:-1] = c[2:] * alfa_i[1:-1] + betta_i[1:-1]
        c[-1] = 0
        c[0] = c[1]

    return c

def main():
    scheme = input('Схема explicit, implicit:')
    c_init = 1
    c = np.linspace(c_init, c_init, num_steps + 2)
    c_final = []
    for t in range(num_steps):
        c = fick(c, scheme)
        c_final.append(c[-2])
        plt.title("Зависимость концентрации от шага")
        plt.xlabel('steps')
        plt.ylabel('concentration')

        if t % 10 == 0:
            plt.plot(c_final)
        #print(c)
    #print(c_fin)
main()
plt.show()