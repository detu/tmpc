#include <tmpc/ocp/SensitivityTuple.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/Testing.hpp>

#include <iostream>

namespace tmpc :: testing
{
	using Kernel = BlazeKernel<double>;

	TEST(SensitivityTupleTest, testCtor)
	{
		SensitivityTuple<Kernel, 2> s {7u, {3u, 4u}};
		//SensitivityTuple<Kernel, 2> s {2u, 3u, 4u};

		EXPECT_EQ(size(s.value()), 7u);
		EXPECT_EQ(rows(s.sens<0>()), 7u);
		EXPECT_EQ(columns(s.sens<0>()), 3u);
		EXPECT_EQ(rows(s.sens<1>()), 7u);
		EXPECT_EQ(columns(s.sens<1>()), 4u);
	}


	/*
	TEST(SensitivityTupleTest, testValue)
	{
		SensitivityTuple<double> s {1.};

		ASSERT_EQ(s.value(), 1.);
	}


	TEST(SensitivityTupleTest, testSens)
	{
		SensitivityTuple<double, double, double> s {1., 2., 3.};

		ASSERT_EQ(s.sens<0>(), 2.);
		ASSERT_EQ(s.sens<1>(), 3.);
	}


	template<class Ch, class Tr, class Tuple, std::size_t... Is>
	void print_tuple_impl(std::basic_ostream<Ch,Tr>& os,
						  const Tuple & t,
						  std::index_sequence<Is...>)
	{
		((os << (Is == 0? "" : ", ") << std::get<Is>(t)), ...);
	}
	 
	template<class Ch, class Tr, class... Args>
	decltype(auto) operator<<(std::basic_ostream<Ch, Tr>& os,
							  const std::tuple<Args...>& t)
	{
		os << "(";
		print_tuple_impl(os, t, std::index_sequence_for<Args...>{});
		return os << ")";
	}


	template<class Tuple1, class Tuple2, std::size_t... Is>
	decltype(auto) add_tuple_impl(const Tuple1 & t1, const Tuple2 & t2,
						  std::index_sequence<Is...>)
	{
		return std::make_tuple(std::get<Is>(t1) + std::get<Is>(t2) ...);
	}


	TEST(SensitivityTupleTest, testAddTuples)
	{
		std::tuple<std::string, double, double> t1 {"a", 2.0, 3.0};		
		std::tuple<std::string, double, float> t2 {"b", 2.0, -3.0f};		
		
		//std::cout << t << std::endl;
		std::cout << add_tuple_impl(t1, t2, std::index_sequence<0, 2> {}) << std::endl;
		//f(t);
		//std::cout <<
	}
	*/
}