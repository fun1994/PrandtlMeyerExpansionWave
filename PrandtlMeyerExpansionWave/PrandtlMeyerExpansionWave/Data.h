#pragma once
#include <vector>

class Data {
public:
	std::vector<double> xi;
	std::vector<double> eta;
	std::vector<double> x;
	std::vector<std::vector<double>> y;
	std::vector<std::vector<double>> rho;
	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> v;
	std::vector<std::vector<double>> T;
	std::vector<std::vector<double>> p;
	std::vector<std::vector<double>> Ma;
	std::vector<std::vector<double>> F1;
	std::vector<std::vector<double>> F2;
	std::vector<std::vector<double>> F3;
	std::vector<std::vector<double>> F4;
};
