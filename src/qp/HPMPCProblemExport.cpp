/*
 * HPMPCProblemExport.cpp
 *
 *  Created on: Jun 28, 2016
 *      Author: kotlyar
 */

#include "HPMPCProblemExport.hpp"
#include <vector>

namespace hpmpc_problem_export
{
	using std::endl;

	void CFormatter::printStructure(int N, int const *nx, int const *nu, int const *nb, int const *ng) const
	{
		_os << "typedef struct {" << endl;
		_os << "int N;" << endl;
		_os << "int k_max;" << endl;
		_os << "double mu0;" << endl;
		_os << "double mu_tol;" << endl;
		_os << "int warm_start;" << endl;

		printArrayStructure("int", "nx", N + 1);
		printArrayStructure("int", "nu", N);
		printArrayStructure("int", "nb", N + 1);
		printArrayStructure("int", "ng", N + 1);

		printArrayOfArraysStructure("double", true, "A", N, nx + 1, nx);
		printArrayOfArraysStructure("double", true, "B", N, nx + 1, nu);
		printArrayOfArraysStructure("double", true, "b", N, nx + 1);
		printArrayOfArraysStructure("double", true, "Q", N + 1, nx, nx);
		printArrayOfArraysStructure("double", true, "S", N, nu, nx);
		printArrayOfArraysStructure("double", true, "R", N, nu, nu);
		printArrayOfArraysStructure("double", true, "q", N + 1, nx);
		printArrayOfArraysStructure("double", true, "r", N, nu);

		{
			std::vector<int> nb(N + 1);
			std::transform(nu, nu + N, nx, nb.begin(), std::plus<int>());
			nb[N] = nx[N];
			printArrayOfArraysStructure("double", true, "lb", N + 1, nb.data());
			printArrayOfArraysStructure("double", true, "ub", N + 1, nb.data());
		}

		printArrayOfArraysStructure("double", true, "C" , N + 1, nx, ng);
		printArrayOfArraysStructure("double", true, "D" , N    , nu, ng);
		printArrayOfArraysStructure("double", true, "lg", N + 1, ng);
		printArrayOfArraysStructure("double", true, "ug", N + 1, ng);

		printArrayOfArraysStructure("double", false, "x", N + 1, nx);
		printArrayOfArraysStructure("double", false, "u", N    , nu);

		_os << "} ProblemStruct;" << endl;
	}

	void CFormatter::printArrayStructure(std::string const& type, std::string const& name, int N) const
	{
		_os << type << " " << name << "[" << N << "];" << endl;
	}

	void CFormatter::printArrayStructure(std::string const& type, std::string const& name, int M, int N) const
	{
		printArrayStructure(type, name, M * N);
	}

	void CFormatter::printArrayOfArraysStructure(std::string const& type, bool const_ptr, std::string const& name, int N, int const * n) const
	{
		printArrayStructure(type + (const_ptr ? " const *" : " *"), name, N);

		for (int i = 0; i < N; ++i)
			printArrayStructure(type, name + std::to_string(i), n[i]);
	}

	void CFormatter::printArrayOfArraysStructure(std::string const& type, bool const_ptr, std::string const& name, int N, int const * m, int const * n) const
	{
		printArrayStructure(type + (const_ptr ? " const *" : " *"), name, N);

		for (int i = 0; i < N; ++i)
			printArrayStructure(type, name + std::to_string(i), m[i], n[i]);
	}

	void CFormatter::printFunctionHeader() const
	{
		_os << "void init_problem(ProblemStruct * " << _varName << ") {" << endl;
	}

	void CFormatter::printFunctionFooter() const
	{
		_os << "}" << endl;
	}
}
