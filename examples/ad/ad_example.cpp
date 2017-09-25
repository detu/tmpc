#include <tmpc/ad/ad.hpp>


#include <json.hpp>

#include <iostream>
#include <tuple>

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

int main(int, char **)
{
    /*
    M = 1    # mass of the cart [kg]
    m = 0.1  # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8  # length of the rod [m]

    p = SX.sym('p')         # horizontal displacement [m]
    theta = SX.sym('theta') # angle with the vertical [rad]
    v = SX.sym('v')         # horizontal velocity [m/s]
    omega = SX.sym('omega') # angular velocity [rad/s]
    F = SX.sym('F')         # horizontal force [N]

    ode_rhs = vertcat(v,
                      omega,
                      (- l*m*sin(theta)*omega**2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)**2),
                      (- l*m*cos(theta)*sin(theta)*omega**2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)**2)))
    */

    Real const x[4] = {0.1, 0.2, 0.3, 0.4};
    Real const u[1] = {0.5};
    Real dx[4];

    std::cout << "x: " << std::endl;
    std::copy_n(x, 4, std::ostream_iterator<Real>(std::cout, ", "));
    std::cout << std::endl;

    std::cout << "u: ";
    std::copy_n(u, 1, std::ostream_iterator<Real>(std::cout, ", "));
    std::cout << std::endl;

    f1(x, u, dx);
    std::cout << "f1(x, u): ";
    std::copy_n(dx, 4, std::ostream_iterator<Real>(std::cout, ", "));
    std::cout << std::endl;

    f2(x, u, dx);
    std::cout << "f2(x, u): ";
    std::copy_n(dx, 4, std::ostream_iterator<Real>(std::cout, ", "));
    std::cout << std::endl;    
    
    return 0;
}