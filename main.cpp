#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

const double alpha = 0.05; // воронка
const double gamma = 0.1;
//const double alpha = 1.1; // аттрактор обычный
//const double gamma = 0.87;
//const double alpha = 0.14; // аттрактор "двойной"
//const double gamma = 0.1;

//const double h = 0.01; // Аттрактор из вики
const double h = 0.001; // "Воронка" из вики
const double T = 100.0;

void f_RF(double * X, double* dX) {
    dX[0] = X[1] * (X[2] - 1 + X[0]*X[0]) + gamma * X[0];
    dX[1] = X[0] * (3 * X[2] + 1 - X[0]*X[0]) + gamma * X[1];
    dX[2] = -2 * X[2] * (alpha + X[0] * X[1]);
}

void three_eights_method(double* state, double* buf, double* k1, double* k2, double* k3, double* k4) {
    f_RF(state, k1);
    for (int i = 0; i < 3; i++) buf[i] = state[i] + h * (k1[i] * 1.0/3);
    f_RF(buf, k2);
    for (int i = 0; i < 3; i++) buf[i] = state[i] + h * (-1.0/3 * k1[i] + 1 * k2[i]);
    f_RF(buf, k3);
    for (int i = 0; i < 3; i++) buf[i] = state[i] + h * (k1[i] - k2[i] + k3[i]);
    f_RF(buf, k4);
    for (int i = 0; i < 3; i++)
        state[i] += h * (1.0 / 8.0 * k1[i] + 3.0 / 8.0 * k2[i] + 3.0 / 8.0 * k3[i] + 1.0 / 8.0 * k4[i]);
}

void Runge_Kutta_4(double* state, double* buf, double* k1, double* k2, double* k3, double* k4) {
    f_RF(state, k1);
    for (int i = 0; i < 3; i++) buf[i] = state[i] + h * 0.5 * k1[i];
    f_RF(buf, k2);
    for (int i = 0; i < 3; i++) buf[i] = state[i] + h * 0.5 * k2[i];
    f_RF(buf, k3);
    for (int i = 0; i < 3; i++) buf[i] = state[i] + h * 1 * k3[i];
    f_RF(buf, k4);
    for (int i = 0; i < 3; i++)
        state[i] += h * (1.0 / 6.0 * k1[i] + 1.0 / 3.0 * k2[i] + 1.0 / 3.0 * k3[i] + 1.0 / 6.0 * k4[i]);
}



int main() {
    double State[3] {0.1, -0.1, 0.1}; // воронка
//    double State[3] {-1.0, 0.0, 0.5}; // аттрактор
    double Buf[3] {};
    double k1[3] {};
    double k2[3] {};
    double k3[3] {};
    double k4[3] {};

//    auto Buf = new double [3];
//    auto k1 = new double [3];
//    auto k2 = new double [3];
//    auto k3 = new double [3];
//    auto k4 = new double [3];

    std::ofstream outputFile("../data.txt");
    for (double t = 0; t < T; t += h) {
        outputFile << std::fixed << std::setprecision(18) << State[0] << " " << State[1]  << " " << State[2] << std::endl;
//        three_eights_method(State, Buf, k1, k2, k3, k4);
        Runge_Kutta_4(State, Buf, k1, k2, k3, k4);
    }
}
