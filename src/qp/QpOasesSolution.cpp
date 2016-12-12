/*
 * QpOasesSolution.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/QpOasesSolution.hpp>

namespace tmpc {

QpOasesSolution::QpOasesSolution(std::vector<QpSize> const& sz)
:	primalSolution_(numVariables(sz))
,	dualSolution_(numEqualities(sz) + numInequalities(sz))
,	size_(sz)
{
	assert(stage_.empty());
	stage_.reserve(sz.size());

	double * x = primalSolution_.data();
	double * lambda_x = dualSolution_.data();
	double * lambda = dualSolution_.data() + primalSolution_.size();

	for (auto s = sz.cbegin(); s != sz.cend(); ++s)
	{
		auto const s_next = s + 1;
		auto const nx_next = s_next != sz.end() ? s_next->nx() : 0;
		stage_.emplace_back(*s, nx_next, x, lambda_x, lambda);

		x += s->nx() + s->nu();
		lambda_x += s->nx() + s->nu();
		lambda += s->nc() + nx_next;
	}
}

QpOasesSolution::Stage::Stage(QpSize const& sz, std::size_t nx_plus, double * x, double * lam_x, double * lam)
:	x_(x, sz.nx())
,	u_(x + sz.nx(), sz.nu())
,	lamX_(lam_x, sz.nx())
,	lamU_(lam_x + sz.nx(), sz.nu())
,	lam_(lam, sz.nc())
,	pi_(lam + sz.nc(), nx_plus)
{
}

}	// namespace tmpc
