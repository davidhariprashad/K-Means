#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <limits>

constexpr double MIN = -100.0;
constexpr double MAX = 100.0;

std::vector<std::pair<double, double>> read_file(const char * filename);
std::vector<size_t> k_means(std::vector<std::pair<double, double>> & points, const int k);
double distance(std::pair<double, double>& a, std::pair<double, double>& b);
template <typename T> std::pair<T, T>& operator+=(std::pair<T, T>& lhs, const std::pair<T, T>& rhs);
template <typename T, typename U> std::pair<T, T>& operator/=(std::pair<T, T>& lhs, const U rhs);
void write_file(const char* filename, const std::vector<std::pair<double, double>>& points, const std::vector<size_t> points_clusters);

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		std::cerr << "<input filename> <k> <output filename>\n";
		return 1;
	}

	int k = std::stoi(argv[2]);
	assert(k > 0);

	auto points = read_file(argv[1]);
	auto points_clusters = k_means(points, k);
	write_file(argv[3], points, points_clusters);

	return 0;
}

std::vector<std::pair<double, double>> read_file(const char * input_filename)
{
	std::vector<std::pair<double, double>> points;
	std::ifstream ifs(input_filename);
	double x, y; char comma;
	while (ifs >> x >> comma >> y)
	{
		points.push_back(std::pair<double, double>(x, y));
	}
	return points;
}

std::vector<size_t> k_means(
	std::vector<std::pair<double, double>> & points,
	const int k)
{
	typedef std::pair<double, double> Point;

	// randomly assign centroids value
	std::vector<Point> centroids(k);
	{
		std::mt19937 engine(std::random_device{}());
		std::uniform_real_distribution<double> distribution(MIN, MAX);
		for (auto& centroid : centroids)
			centroid = Point(distribution(engine), distribution(engine));
	}

	// assign cluster for each point
	std::vector<size_t> points_clusters(points.size());
	std::vector<size_t> clusters_sizes(k);
	std::fill(clusters_sizes.begin(), clusters_sizes.end(), 0U);
	for (size_t point = 0U; point < points.size(); ++point)
	{
		double minimum_distance = distance(points[point], centroids[0]);
		size_t closest_centroid = 0U;
		for (size_t centroid = 1U; centroid < centroids.size(); ++centroid)
		{
			const double d = distance(points[point], centroids[centroid]);
			if (d < minimum_distance)
			{
				minimum_distance = d;
				closest_centroid = centroid;
			}
		}
		points_clusters[point] = closest_centroid;
		++clusters_sizes[closest_centroid];
	}

	// do while change, recalculate centroids and reassign points to clusters
	int changes;
	do
	{
		changes = 0;

		// recalculate centroids
		for (auto& centroid : centroids)
			centroid = {0.0, 0.0};

		for (size_t point = 0U; point < points.size(); ++point)
			centroids[points_clusters[point]] += {points[point].first, points[point].second};

		for (size_t centroid = 0U; centroid < centroids.size(); ++centroid)
			if (clusters_sizes[centroid] != 0U)
				centroids[centroid] /= clusters_sizes[centroid];

		// reassign points to clusters
		for (size_t point = 0U; point < points.size(); ++point)
		{
			double minimum_distance = distance(points[point], centroids[0]);
			size_t closest_centroid = 0U;
			for (size_t centroid = 1U; centroid < centroids.size(); ++centroid)
			{
				double d = distance(points[point], centroids[centroid]);
				if (d < minimum_distance)
				{
					minimum_distance = d;
					closest_centroid = centroid;
				}
			}
			if (points_clusters[point] != closest_centroid)
			{
				++changes;
				++clusters_sizes[closest_centroid];
				--clusters_sizes[points_clusters[point]];
				points_clusters[point] = closest_centroid;
			}
		}

		std::cout << "changes: " << changes << '\n';
	} while (changes > 0);

	std::cout << "Centroids:\n";
	for (auto& centroid : centroids)
		std::cout << centroid.first << ',' << centroid.second << '\n';

	return points_clusters;
}

double distance(
	std::pair<double, double>& a,
	std::pair<double, double>& b)
{
	return std::sqrt(std::pow(a.first - b.first, 2.0) + std::pow(a.second - b.second, 2));
}

template <typename T>
std::pair<T, T> & operator+=(std::pair<T, T> & lhs, const std::pair<T, T> & rhs)
{
	lhs.first += rhs.first;
	lhs.second += rhs.second;
	return lhs;
}

template <typename T, typename U>
std::pair<T, T> & operator/=(std::pair<T, T>& lhs, const U rhs)
{
	lhs.first /= rhs;
	lhs.second /= rhs;
	return lhs;
}

void write_file(
	const char* filename,
	const std::vector<std::pair<double, double>>& points,
	const std::vector<size_t> points_clusters)
{
	std::ofstream ofs(filename);
	for (size_t point = 0U; point < points.size(); ++point)
	{
		ofs << points[point].first << ',';
		ofs << points[point].second << ',';
		ofs << points_clusters[point] << '\n';
	}
}
