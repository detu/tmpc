#include <tmpc/ad/ad.hpp>


#include <json.hpp>

#include <benchmark/benchmark.h>

#include <iostream>
#include <tuple>
#include <chrono>

using nlohmann::json;

using Real = double;

Real constexpr M = 1.;      // mass of the cart [kg]
Real constexpr m = 0.1;     // mass of the ball [kg]
Real constexpr g = 9.81;    // gravity constant [m/s^2]
Real constexpr l = 0.8;     // length of the rod [m]

namespace tmpc
{
    template <typename OP, typename L, typename R>
    void to_json(json& j, const BinaryOp<OP, L, R>& expr) 
    {
        j = json {{OP::name(), {expr.left(), expr.right()}}};
    }


    template <typename OP, typename A>
    void to_json(json& j, const UnaryOp<OP, A>& expr) 
    {
        j = json {{OP::name(), expr.arg()}};
    }


    template <std::size_t N>
    void to_json(json& j, const Variable<N>& var) 
    {
        j = "@" + std::to_string(N);
    }
}

void f1(Real const x[4], Real const u[1], Real dx[4])
{
    auto const p = x[0];
    auto const theta = x[1];
    auto const v = x[2];
    auto const omega = x[3];
    auto const F = u[0];

    dx[0] = v;
    dx[1] = omega;
    dx[2] = (-l * m * sin(theta) * pow(omega, 2) + F + g * m * cos(theta) * sin(theta)) / (M + m - m * pow(cos(theta), 2));
    dx[3] = (-l * m * cos(theta) * sin(theta) * pow(omega, 2) + F * cos(theta) + g * m * sin(theta) + M * g * sin(theta)) 
        / (l * (M + m - m * pow(cos(theta), 2)));
}

void f2(Real const x[4], Real const u[1], Real dx[4])
{
    using tmpc::Variable;

    Variable<0> p;      // horizontal displacement [m]
    Variable<1> theta;  // angle with the vertical [rad]
    Variable<2> v;      // horizontal velocity [m/s]
    Variable<3> omega;  // angular velocity [rad/s]
    Variable<4> F;      // horizontal force [N]
    
    auto p_dot = v;
    auto theta_dot = omega;
    auto v_dot = (-l * m * sin(theta) * pow(omega, 2) + F + g * m * cos(theta) * sin(theta)) / (M + m - m * pow(cos(theta), 2));
    auto omega_dot = (-l * m * cos(theta) * sin(theta) * pow(omega, 2) + F * cos(theta) + g * m * sin(theta) + M * g * sin(theta)) 
        / (l * (M + m - m * pow(cos(theta), 2)));

    auto arg = std::make_tuple(x[0], x[1], x[2], x[3], u[0]);

    dx[0] = eval(p_dot, arg);
    dx[1] = eval(theta_dot, arg);
    dx[2] = eval(v_dot, arg);
    dx[3] = eval(omega_dot, arg);
}


static void BM_f1(benchmark::State& state) 
{
    Real const x[4] = {0.1, 0.2, 0.3, 0.4};
    Real const u[1] = {0.5};
    Real dx[4];

    while (state.KeepRunning())
        f1(x, u, dx);
}

BENCHMARK(BM_f1);


static void BM_f2(benchmark::State& state) 
{
    Real const x[4] = {0.1, 0.2, 0.3, 0.4};
    Real const u[1] = {0.5};
    Real dx[4];

    while (state.KeepRunning())
        f2(x, u, dx);
}

BENCHMARK(BM_f2);

BENCHMARK_MAIN();