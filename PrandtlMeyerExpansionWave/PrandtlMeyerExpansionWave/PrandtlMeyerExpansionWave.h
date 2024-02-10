#pragma once
#include <iostream>
#include "Data.h"
#define PI 3.141592653589793

class PrandtlMeyerExpansionWave {
	double gamma;
	double R;
	double L;
	double H;
	double E;
	double theta;
	double rho0;
	double T0;
	double Ma0;
	int N;
	double deta;
	double C;
	double Cy;
	double tol;
	void initialize(Data& data, std::vector<double>& y, std::vector<double>& rho, std::vector<double>& u, std::vector<double>& v, std::vector<double>& T, std::vector<double>& p, std::vector<double>& Ma, std::vector<double>& F1, std::vector<double>& F2, std::vector<double>& F3, std::vector<double>& F4) {
		data.eta.resize(N + 1);
		for (int i = 0; i < N + 1; i++) {
			data.eta[i] = i * deta;
		}
		for (int i = 0; i < N + 1; i++) {
			y[i] = data.eta[i] * H;
		}
		for (int i = 0; i < N + 1; i++) {
			rho[i] = rho0;
		}
		for (int i = 0; i < N + 1; i++) {
			u[i] = Ma0 * sqrt(gamma * R * T0);
		}
		for (int i = 0; i < N + 1; i++) {
			v[i] = 0;
		}
		for (int i = 0; i < N + 1; i++) {
			T[i] = T0;
		}
		for (int i = 0; i < N + 1; i++) {
			p[i] = rho0 * R * T0;
		}
		for (int i = 0; i < N + 1; i++) {
			Ma[i] = Ma0;
		}
		for (int i = 0; i < N + 1; i++) {
			F1[i] = rho[i] * u[i];
		}
		for (int i = 0; i < N + 1; i++) {
			F2[i] = rho[i] * pow(u[i], 2) + p[i];
		}
		for (int i = 0; i < N + 1; i++) {
			F3[i] = rho[i] * u[i] * v[i];
		}
		for (int i = 0; i < N + 1; i++) {
			F4[i] = gamma / (gamma - 1) * p[i] * u[i] + rho[i] * u[i] * (pow(u[i], 2) + pow(v[i], 2)) / 2;
		}
		data.xi.push_back(0);
		data.x.push_back(0);
		data.y.push_back(y);
		data.rho.push_back(rho);
		data.u.push_back(u);
		data.v.push_back(v);
		data.T.push_back(T);
		data.p.push_back(p);
		data.Ma.push_back(Ma);
		data.F1.push_back(F1);
		data.F2.push_back(F2);
		data.F3.push_back(F3);
		data.F4.push_back(F4);
	}
	double stepSize(double h, std::vector<double>& Ma) {
		std::vector<double> mu(N + 1), tan1(N + 1), tan2(N + 1);
		for (int i = 0; i < N + 1; i++) {
			mu[i] = asin(1 / Ma[i]);
		}
		for (int i = 0; i < N + 1; i++) {
			tan1[i] = abs(tan(theta + mu[i]));
		}
		for (int i = 0; i < N + 1; i++) {
			tan2[i] = abs(tan(theta - mu[i]));
		}
		return C * deta * h / std::max(*std::max_element(tan1.begin(), tan1.end()), *std::max_element(tan2.begin(), tan2.end()));
	}
	double f(double Ma) {
		return sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (pow(Ma, 2) - 1))) - atan(sqrt(pow(Ma, 2) - 1));
	}
	double dfdMa(double Ma) {
		return (Ma / (1 + (gamma - 1) / (gamma + 1) * (pow(Ma, 2) - 1)) - 1 / Ma) / sqrt(pow(Ma, 2) - 1);
	}
	double Newton(double f0, double Ma0) {
		double Ma = Ma0;
		while (abs(f(Ma) - f0) > tol) {
			Ma -= (f(Ma) - f0) / dfdMa(Ma);
		}
		return Ma;
	}
public:
	PrandtlMeyerExpansionWave(double gamma, double R0, double M, double L, double H, double E, double theta, double rho0, double T0, double Ma0, int N, double C, double Cy, double tol) :gamma(gamma), R(R0 / M), L(L), H(H), E(E), theta(theta / 180 * PI), rho0(rho0), T0(T0), Ma0(Ma0), N(N), deta(1.0 / N), C(C), Cy(Cy), tol(tol) {}
	void MacCormack(Data& data) {
		double h = H;
		std::vector<double> y(N + 1), rho(N + 1), u(N + 1), v(N + 1), T(N + 1), p(N + 1), Ma(N + 1), F1(N + 1), F2(N + 1), F3(N + 1), F4(N + 1), detadx(N + 1), G1(N + 1), G2(N + 1), G3(N + 1), G4(N + 1), dF1dxi(N + 1), dF2dxi(N + 1), dF3dxi(N + 1), dF4dxi(N + 1), F1_bar(N + 1), F2_bar(N + 1), F3_bar(N + 1), F4_bar(N + 1), A(N + 1), B(N + 1), C(N + 1), rho_bar(N + 1), p_bar(N + 1), G1_bar(N + 1), G2_bar(N + 1), G3_bar(N + 1), G4_bar(N + 1), dF1dxi_bar(N + 1), dF2dxi_bar(N + 1), dF3dxi_bar(N + 1), dF4dxi_bar(N + 1), dF1dxi_av(N + 1), dF2dxi_av(N + 1), dF3dxi_av(N + 1), dF4dxi_av(N + 1), SF1(N + 1), SF2(N + 1), SF3(N + 1), SF4(N + 1), SF1_bar(N + 1), SF2_bar(N + 1), SF3_bar(N + 1), SF4_bar(N + 1);
		initialize(data, y, rho, u, v, T, p, Ma, F1, F2, F3, F4);
		do {
			double dxi = stepSize(h, Ma);
			for (int i = 0; i < N + 1; i++) {
				detadx[i] = data.xi[data.xi.size() - 1] < E ? 0 : (1 - data.eta[i]) / h * tan(theta);
			}
			for (int i = 0; i < N + 1; i++) {
				G1[i] = rho[i] * F3[i] / F1[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G2[i] = F3[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G3[i] = rho[i] * pow(F3[i] / F1[i], 2) + F2[i] - pow(F1[i], 2) / rho[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G4[i] = gamma / (gamma - 1) * (F2[i] - pow(F1[i], 2) / rho[i]) * F3[i] / F1[i] + rho[i] / 2 * F3[i] / F1[i] * (pow(F1[i] / rho[i], 2) + pow(F3[i] / F1[i], 2));
			}
			for (int i = 0; i < N; i++) {
				dF1dxi[i] = -detadx[i] * (F1[i + 1] - F1[i]) / deta - 1 / h * (G1[i + 1] - G1[i]) / deta;
			}
			for (int i = 0; i < N; i++) {
				dF2dxi[i] = -detadx[i] * (F2[i + 1] - F2[i]) / deta - 1 / h * (G2[i + 1] - G2[i]) / deta;
			}
			for (int i = 0; i < N; i++) {
				dF3dxi[i] = -detadx[i] * (F3[i + 1] - F3[i]) / deta - 1 / h * (G3[i + 1] - G3[i]) / deta;
			}
			for (int i = 0; i < N; i++) {
				dF4dxi[i] = -detadx[i] * (F4[i + 1] - F4[i]) / deta - 1 / h * (G4[i + 1] - G4[i]) / deta;
			}
			for (int i = 1; i < N; i++) {
				SF1[i] = Cy * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (F1[i + 1] - 2 * F1[i] + F1[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				SF2[i] = Cy * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (F2[i + 1] - 2 * F2[i] + F2[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				SF3[i] = Cy * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (F3[i + 1] - 2 * F3[i] + F3[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				SF4[i] = Cy * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (F4[i + 1] - 2 * F4[i] + F4[i - 1]);
			}
			for (int i = 0; i < N; i++) {
				F1_bar[i] = F1[i] + dF1dxi[i] * dxi + SF1[i];
			}
			F1_bar[N] = data.F1[0][N];
			for (int i = 0; i < N; i++) {
				F2_bar[i] = F2[i] + dF2dxi[i] * dxi + SF2[i];
			}
			F2_bar[N] = data.F2[0][N];
			for (int i = 0; i < N; i++) {
				F3_bar[i] = F3[i] + dF3dxi[i] * dxi + SF3[i];
			}
			F3_bar[N] = data.F3[0][N];
			for (int i = 0; i < N; i++) {
				F4_bar[i] = F4[i] + dF4dxi[i] * dxi + SF4[i];
			}
			F4_bar[N] = data.F4[0][N];
			for (int i = 0; i < N + 1; i++) {
				A[i] = pow(F3_bar[i], 2) / (2 * F1_bar[i]) - F4_bar[i];
			}
			for (int i = 0; i < N + 1; i++) {
				B[i] = gamma / (gamma - 1) * F1_bar[i] * F2_bar[i];
			}
			for (int i = 0; i < N + 1; i++) {
				C[i] = -(gamma + 1) / (2 * (gamma - 1)) * pow(F1_bar[i], 3);
			}
			for (int i = 0; i < N + 1; i++) {
				rho_bar[i] = (-B[i] + sqrt(pow(B[i], 2) - 4 * A[i] * C[i])) / (2 * A[i]);
			}
			for (int i = 0; i < N + 1; i++) {
				p_bar[i] = F2_bar[i] - pow(F1_bar[i], 2) / rho_bar[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G1_bar[i] = rho_bar[i] * F3_bar[i] / F1_bar[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G2_bar[i] = F3_bar[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G3_bar[i] = rho_bar[i] * pow(F3_bar[i] / F1_bar[i], 2) + F2_bar[i] - pow(F1_bar[i], 2) / rho_bar[i];
			}
			for (int i = 0; i < N + 1; i++) {
				G4_bar[i] = gamma / (gamma - 1) * (F2_bar[i] - pow(F1_bar[i], 2) / rho_bar[i]) * F3_bar[i] / F1_bar[i] + rho_bar[i] / 2 * F3_bar[i] / F1_bar[i] * (pow(F1_bar[i] / rho_bar[i], 2) + pow(F3_bar[i] / F1_bar[i], 2));
			}
			dF1dxi_bar[0] = -detadx[0] * (F1_bar[1] - F1_bar[0]) / deta - 1 / h * (G1_bar[1] - G1_bar[0]) / deta;
			for (int i = 1; i < N + 1; i++) {
				dF1dxi_bar[i] = -detadx[i] * (F1_bar[i] - F1_bar[i - 1]) / deta - 1 / h * (G1_bar[i] - G1_bar[i - 1]) / deta;
			}
			dF2dxi_bar[0] = -detadx[0] * (F2_bar[1] - F2_bar[0]) / deta - 1 / h * (G2_bar[1] - G2_bar[0]) / deta;
			for (int i = 1; i < N + 1; i++) {
				dF2dxi_bar[i] = -detadx[i] * (F2_bar[i] - F2_bar[i - 1]) / deta - 1 / h * (G2_bar[i] - G2_bar[i - 1]) / deta;
			}
			dF3dxi_bar[0] = -detadx[0] * (F3_bar[1] - F3_bar[0]) / deta - 1 / h * (G3_bar[1] - G3_bar[0]) / deta;
			for (int i = 1; i < N + 1; i++) {
				dF3dxi_bar[i] = -detadx[i] * (F3_bar[i] - F3_bar[i - 1]) / deta - 1 / h * (G3_bar[i] - G3_bar[i - 1]) / deta;
			}
			dF4dxi_bar[0] = -detadx[0] * (F4_bar[1] - F4_bar[0]) / deta - 1 / h * (G4_bar[1] - G4_bar[0]) / deta;
			for (int i = 1; i < N + 1; i++) {
				dF4dxi_bar[i] = -detadx[i] * (F4_bar[i] - F4_bar[i - 1]) / deta - 1 / h * (G4_bar[i] - G4_bar[i - 1]) / deta;
			}
			for (int i = 0; i < N; i++) {
				dF1dxi_av[i] = (dF1dxi[i] + dF1dxi_bar[i]) / 2;
			}
			for (int i = 0; i < N; i++) {
				dF2dxi_av[i] = (dF2dxi[i] + dF2dxi_bar[i]) / 2;
			}
			for (int i = 0; i < N; i++) {
				dF3dxi_av[i] = (dF3dxi[i] + dF3dxi_bar[i]) / 2;
			}
			for (int i = 0; i < N; i++) {
				dF4dxi_av[i] = (dF4dxi[i] + dF4dxi_bar[i]) / 2;
			}
			for (int i = 1; i < N; i++) {
				SF1_bar[i] = Cy * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (F1_bar[i + 1] - 2 * F1_bar[i] + F1_bar[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				SF2_bar[i] = Cy * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (F2_bar[i + 1] - 2 * F2_bar[i] + F2_bar[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				SF3_bar[i] = Cy * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (F3_bar[i + 1] - 2 * F3_bar[i] + F3_bar[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				SF4_bar[i] = Cy * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (F4_bar[i + 1] - 2 * F4_bar[i] + F4_bar[i - 1]);
			}
			for (int i = 0; i < N; i++) {
				F1[i] += dF1dxi_av[i] * dxi + SF1_bar[i];
			}
			for (int i = 0; i < N; i++) {
				F2[i] += dF2dxi_av[i] * dxi + SF2_bar[i];
			}
			for (int i = 0; i < N; i++) {
				F3[i] += dF3dxi_av[i] * dxi + SF3_bar[i];
			}
			for (int i = 0; i < N; i++) {
				F4[i] += dF4dxi_av[i] * dxi + SF4_bar[i];
			}
			F1[N] = data.F1[0][N];
			F2[N] = data.F2[0][N];
			F3[N] = data.F3[0][N];
			F4[N] = data.F4[0][N];
			for (int i = 0; i < N + 1; i++) {
				A[i] = pow(F3[i], 2) / (2 * F1[i]) - F4[i];
			}
			for (int i = 0; i < N + 1; i++) {
				B[i] = gamma / (gamma - 1) * F1[i] * F2[i];
			}
			for (int i = 0; i < N + 1; i++) {
				C[i] = -(gamma + 1) / (2 * (gamma - 1)) * pow(F1[i], 3);
			}
			for (int i = 0; i < N + 1; i++) {
				rho[i] = (-B[i] + sqrt(pow(B[i], 2) - 4 * A[i] * C[i])) / (2 * A[i]);
			}
			for (int i = 0; i < N + 1; i++) {
				u[i] = F1[i] / rho[i];
			}
			for (int i = 0; i < N + 1; i++) {
				v[i] = F3[i] / F1[i];
			}
			for (int i = 0; i < N + 1; i++) {
				p[i] = F2[i] - F1[i] * u[i];
			}
			for (int i = 0; i < N + 1; i++) {
				T[i] = p[i] / (rho[i] * R);
			}
			for (int i = 0; i < N + 1; i++) {
				Ma[i] = sqrt(pow(u[i], 2) + pow(v[i], 2)) / sqrt(gamma * R * T[i]);
			}
			data.xi.push_back(data.xi[data.xi.size() - 1] + dxi);
			double phi = data.xi[data.xi.size() - 1] < E ? atan(v[0] / u[0]) : theta + atan(v[0] / u[0]);
			double Ma_cal = Ma[0];
			double f_cal = f(Ma_cal);
			double f_act = f_cal + phi;
			double Ma_act = Newton(f_act, Ma_cal);
			Ma[0] = Ma_act;
			p[0] *= pow((1 + (gamma - 1) / 2 * pow(Ma_cal, 2)) / (1 + (gamma - 1) / 2 * pow(Ma_act, 2)), gamma / (gamma - 1));
			T[0] *= (1 + (gamma - 1) / 2 * pow(Ma_cal, 2)) / (1 + (gamma - 1) / 2 * pow(Ma_act, 2));
			rho[0] = p[0] / (R * T[0]);
			u[0] = data.xi[data.xi.size() - 1] < E ? Ma[0] * sqrt(gamma * R * T[0]) : Ma[0] * sqrt(gamma * R * T[0]) * cos(theta);
			v[0] = data.xi[data.xi.size() - 1] < E ? 0 : -Ma[0] * sqrt(gamma * R * T[0]) * sin(theta);
			F1[0] = rho[0] * u[0];
			F2[0] = rho[0] * pow(u[0], 2) + p[0];
			F3[0] = rho[0] * u[0] * v[0];
			F4[0] = gamma / (gamma - 1) * p[0] * u[0] + rho[0] * u[0] * (pow(u[0], 2) + pow(v[0], 2)) / 2;
			data.x.push_back(data.xi[data.xi.size() - 1]);
			h = data.xi[data.xi.size() - 1] < E ? H : H + (data.xi[data.xi.size() - 1] - E) * tan(theta);
			for (int i = 0; i < N + 1; i++) {
				y[i] = data.xi[data.xi.size() - 1] < E ? data.eta[i] * h : data.eta[i] * h - (data.xi[data.xi.size() - 1] - E) * tan(theta);
			}
			data.y.push_back(y);
			data.rho.push_back(rho);
			data.u.push_back(u);
			data.v.push_back(v);
			data.T.push_back(T);
			data.p.push_back(p);
			data.Ma.push_back(Ma);
			data.F1.push_back(F1);
			data.F2.push_back(F2);
			data.F3.push_back(F3);
			data.F4.push_back(F4);
		} while (data.xi[data.xi.size() - 1] < L);
	}
};
