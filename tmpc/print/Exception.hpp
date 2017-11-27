#pragma once

#include <stdexcept>
#include <ostream>


namespace tmpc
{
    // prints the explanatory string of an exception. If the exception is nested,
    // recurses to print the explanatory of the exception it holds
    void print(std::ostream& os, const std::exception& e, int level = 0);
}