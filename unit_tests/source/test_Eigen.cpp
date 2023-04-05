#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>


#include "simple_math/Utility.h"
#include "simple_math/EigenDefines.h"
#include "simple_math/HeapPreallocatedMatrix.h"

using namespace simple_math;

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
	template<typename F, typename ...Args>
	static typename TimeT::rep execution(F&& func, Args&&... args)
	{
		auto start = std::chrono::steady_clock::now();
		std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
		auto duration = std::chrono::duration_cast<TimeT>
			(std::chrono::steady_clock::now() - start);
		return duration.count();
	}
};

void DynamicResize(int n)
{
	GeneralMatrix mat;
	for (auto i = 0; i < n; i++)
	{
		mat.resize(50, 50);
		mat.setOnes();
		mat.resize(120, 120);
		mat.setOnes();
	}
}

void HeapPreallocatedResize(int n)
{
	HeapPreallocatedMatrix<double, 120, 120> mat;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 120, 120> mat2;
	for (auto i = 0; i < n; i++)
	{
		mat.resize(50, 50);
		mat.setOnes();
		mat.resize(120, 120);
		mat.setOnes();
	}
}

void HeapPreallocatedResizeWithAssignment(int n)
{
	HeapPreallocatedMatrix<double, 120, 120> mat;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 120, 120> mat2;
	for (auto i = 0; i < n; i++)
	{
		mat.resize(50, 50);
		mat.setOnes();

		mat2.resize(50, 50);
		mat2 = mat;

		mat.resize(120, 120);
		mat.setOnes();

		mat2.resize(120, 120);
		mat2 = mat;
	}
}

void StackPreallocatedResize(int n)
{
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 120, 120> mat;
	for (auto i = 0; i < n; i++)
	{
		mat.resize(50, 50);
		mat.setOnes();
		mat.resize(120, 120);
		mat.setOnes();
	}
}

int main()
{
	int n = 1000;

	auto bench_dyn = [n]() {
		DynamicResize(n);
	};
	std::cout << "Time for resizing dynamically reallocated matrices: "
		<< measure<>::execution(bench_dyn)
		<< " milliseconds." << std::endl;

	auto bench_prealloc = [n]() {
		HeapPreallocatedResize(n);
	};
	std::cout << "Time for resizing heap preallocated matrices: "
		<< measure<>::execution(bench_prealloc)
		<< " milliseconds." << std::endl;

	auto bench_prealloc2 = [n]() {
		StackPreallocatedResize(n);
	};
	std::cout << "Time for resizing stack preallocated matrices: "
		<< measure<>::execution(bench_prealloc2)
		<< " milliseconds." << std::endl;

	
	return 0;
}