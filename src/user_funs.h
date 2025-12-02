#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);

matrix ff3T(matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);
matrix df3(double, matrix, matrix = NAN, matrix = NAN);

matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix Hf4T(matrix, matrix = NAN, matrix = NAN);

// Lab 4 - Klasyfikator logistyczny
double sigmoid(double);
double h_theta(matrix, matrix);
matrix ff4R(matrix, matrix = NAN, matrix = NAN);
matrix gf4R(matrix, matrix = NAN, matrix = NAN);

