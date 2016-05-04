#pragma once

#include <stdexcept>

#include <qpOASES.hpp>

namespace qpOASES
{
	struct Exception : std::runtime_error
	{
		Exception(returnValue ret)
		: std::runtime_error(ErrorMessage(ret)), _returnValue(ret)
		{}

		returnValue getReturnValue() const { return _returnValue; }

	private:
		returnValue _returnValue;

		static std::string ErrorMessage(returnValue ret)
		{
			return "qpOASES return code " + std::to_string(ret);
		}
	};
}
