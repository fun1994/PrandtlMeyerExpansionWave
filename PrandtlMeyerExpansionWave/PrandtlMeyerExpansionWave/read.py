# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 11:04:52 2024

@author: acer
"""

import numpy as np


def read_1d(filename):
    with open("./data/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read_2d(filename):
    data = []
    with open("./data/" + filename + ".txt", "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            data_temp = line.split()
            data.append(data_temp)
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = float(data[i][j])
    data = np.array(data)
    return data

def read():
    xi = read_1d("xi")
    eta = read_1d("eta")
    x = read_1d("x")
    y = read_2d("y")
    rho = read_2d("rho")
    u = read_2d("u")
    v = read_2d("v")
    T = read_2d("T")
    p = read_2d("p")
    Ma = read_2d("Ma")
    F1 = read_2d("F1")
    F2 = read_2d("F2")
    F3 = read_2d("F3")
    F4 = read_2d("F4")
    return xi, eta, x, y, rho, u, v, T, p, Ma, F1, F2, F3, F4


xi, eta, x, y, rho, u, v, T, p, Ma, F1, F2, F3, F4 = read()
