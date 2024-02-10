#include <fstream>
#include "PrandtlMeyerExpansionWave.h"

void save(std::vector<double>& data, std::string filename) {
	std::ofstream file("./data/" + filename + ".txt");
	for (int i = 0; i < data.size(); i++) {
		file << data[i];
		if (i < data.size() - 1) {
			file << " ";
		}
	}
	file.close();
}

void save(std::vector<std::vector<double>>& data, std::string filename) {
	std::ofstream file("./data/" + filename + ".txt");
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data[i].size(); j++) {
			file << data[i][j];
			if (j < data[i].size() - 1) {
				file << " ";
			}
		}
		if (i < data.size() - 1) {
			file << "\n";
		}
	}
	file.close();
}

void save(Data& data) {
	save(data.xi, "xi");
	save(data.eta, "eta");
	save(data.x, "x");
	save(data.y, "y");
	save(data.rho, "rho");
	save(data.u, "u");
	save(data.v, "v");
	save(data.T, "T");
	save(data.p, "p");
	save(data.Ma, "Ma");
	save(data.F1, "F1");
	save(data.F2, "F2");
	save(data.F3, "F3");
	save(data.F3, "F4");
}

int main() {
	PrandtlMeyerExpansionWave PMEW(1.4, 8.314, 0.029, 65, 40, 10, 5.325, 1.23, 286.1, 2, 40, 0.5, 0.6, 1e-6);
	Data data;
	PMEW.MacCormack(data);
	save(data);
	return 0;
}
