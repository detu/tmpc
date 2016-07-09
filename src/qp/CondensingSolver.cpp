#include <qp/CondensingSolver.hpp>

#include <sstream>
#include <stdexcept>

static qpOASES::Options& ensureConsistency(qpOASES::Options& options)
{
	auto const ret = options.ensureConsistency();
	if(ret != qpOASES::SUCCESSFUL_RETURN)
	{
		std::stringstream msg;
		msg << "ensureConsistency(): qpOASES::Options::ensureConsistency() returned error code " << ret << ".";
		throw std::runtime_error(msg.str());
	}

	return options;
}

static qpOASES::Options& setToReliable(qpOASES::Options& options)
{
	auto const ret = options.setToReliable();
	if(ret != qpOASES::SUCCESSFUL_RETURN)
	{
		std::stringstream msg;
		msg << "setToReliable(): qpOASES::Options::setToReliable() returned error code " << ret << ".";
		throw std::runtime_error(msg.str());
	}

	return options;
}

static qpOASES::Options& setToMPC(qpOASES::Options& options)
{
	auto const ret = options.setToMPC();
	if(ret != qpOASES::SUCCESSFUL_RETURN)
	{
		std::stringstream msg;
		msg << "setToMPC(): qpOASES::Options::setToMPC() returned error code " << ret << ".";
		throw std::runtime_error(msg.str());
	}

	return options;
}

namespace tmpc
{
	qpOASESOptions::qpOASESOptions(qpOASES::Options const& options)
	:	_options()
	{
		ensureConsistency(_options);
	}

	qpOASESOptions::operator qpOASES::Options const&() const
	{
		return _options;
	}

	qpOASESOptions qpOASESOptions::Reliable()
	{
		qpOASES::Options options;
		return setToReliable(options);
	}

	qpOASESOptions qpOASESOptions::MPC()
	{
		qpOASES::Options options;
		return setToMPC(options);
	}

	qpOASES::PrintLevel qpOASESOptions::getPrintLevel() const
	{
		return _options.printLevel;
	}

	qpOASESOptions& qpOASESOptions::setPrintLevel(qpOASES::PrintLevel level)
	{
		_options.printLevel = level;
		ensureConsistency(_options);

		return *this;
	}
}
