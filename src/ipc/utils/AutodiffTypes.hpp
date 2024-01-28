#pragma once

#include "eigen_ext.hpp"
#include "autodiff.h"

namespace ipc
{
	template <int dim, int max_dim = dim>
	using AutodiffScalarGrad = DScalar1<double, Eigen::Matrix<double, dim, 1, 0, max_dim, 1>>;
	template <int dim, int max_dim = dim>
	using AutodiffScalarHessian = DScalar2<double, Eigen::Matrix<double, dim, 1, 0, max_dim, 1>, Eigen::Matrix<double, dim, dim, 0, max_dim, max_dim>>;

	template <class T>
	class AutoDiffAllocator
	{
	public:
		T operator()(const int i, const double &v) const
		{
			return T(i, v);
		}
	};

	template <>
	class AutoDiffAllocator<double>
	{
	public:
		double operator()(const int i, const double &v) const
		{
			return v;
		}
	};
} // namespace polyfem
