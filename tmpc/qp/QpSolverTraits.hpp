#pragma once


namespace tmpc
{
   	/// @brief Determines whether a specific solver supports solving tree problems.
	///
	/// By default it is assumed that it does. If it does not,
	/// specialize the SupportsTreeProblems for your solver class.
	///
	template <typename Solver>
	struct SupportsTreeProblems
	{
		static constexpr bool value = true;
	};
}