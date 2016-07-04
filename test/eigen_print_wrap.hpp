/*
 * Makes the Eigen3 matrix classes printable from GTest checking macros like EXPECT_EQ.
 * Usage: EXPECT_EQ(print_wrap(a), print_wrap(b))
 * Taken from this post: http://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */

#pragma once

#include <iostream>

template <class Base>
class EigenPrintWrap : public Base {
    friend void PrintTo(const EigenPrintWrap &m, ::std::ostream *o) {
        *o << "\n" << m;
    }
};

template <class Base>
const EigenPrintWrap<Base> &print_wrap(const Base &base) {
    return static_cast<const EigenPrintWrap<Base> &>(base);
}
