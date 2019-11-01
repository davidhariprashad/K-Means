#include <fstream>
#include <string>
#include <random>
#include <utility>
#include <cmath>
#include <cassert>

constexpr double PI = 3.14159265358979323846;
constexpr double RADIUS = 25.0;

int main(int argc, char * argv[])
{
	assert(argc == 3);
	const int N = std::stoi(argv[2]);
	assert(N > 1);

	const std::pair<double, double> c1(-25.0, 15.0);
	const std::pair<double, double> c2(20.0, -20.0);

	std::mt19937 engine(std::random_device{}());
	std::uniform_real_distribution<double> distance_distribution(0.0, RADIUS);
	std::uniform_real_distribution<double> radian_distribution(0.0, 2.0*PI);
	std::ofstream ofs(argv[1]);
	for (int p = 0; p < N; ++p)
	{
		double angle = radian_distribution(engine);
		double distance = distance_distribution(engine);
		ofs << c1.first + std::cos(angle)*distance << ',';
		ofs << c1.second + std::sin(angle)*distance << '\n';
	}
	for (int p = 0; p < N; ++p)
	{
		double angle = radian_distribution(engine);
		double distance = distance_distribution(engine);
		ofs << c2.first + std::cos(angle)*distance << ',';
		ofs << c2.second + std::sin(angle)*distance << '\n';
	}
	return 0;
}
