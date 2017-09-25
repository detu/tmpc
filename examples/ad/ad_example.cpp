#include <tmpc/ad/ad.hpp>


#include <json.hpp>

#include <iostream>
#include <tuple>
#include <chrono>

using nlohmann::json;

using Real = double;

Real constexpr M = 1.;      // mass of the cart [kg]
Real constexpr m = 0.1;     // mass of the ball [kg]
Real constexpr g = 9.81;    // gravity constant [m/s^2]
Real constexpr l = 0.8;     // length of the rod [m]

volatile int vv {};
volatile Real result {};

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
        j = json {};
        j[OP::name()] = expr.arg();
    }


    template <std::size_t N>
    void to_json(json& j, const Variable<N>& var) 
    {
        j = "@" + std::to_string(N);
    }


    void to_json(json& j, const Zero& expr) 
    {
        j = "ZERO";
    }


    void to_json(json& j, const Identity& expr) 
    {
        j = "IDENTITY";
    }
}


template <typename T>
std::ostream& operator<<(std::ostream& os, tmpc::ExpressionBase<T> const& expr)
{
    json j = expr.derived();
    return os << j.dump(4);
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

    ++vv;
}

void f2(Real const x[4], Real const u[1], Real dx[4])
{
    using tmpc::Variable;

    Variable<0> constexpr p {};      // horizontal displacement [m]
    Variable<1> constexpr theta {};  // angle with the vertical [rad]
    Variable<2> constexpr v {};      // horizontal velocity [m/s]
    Variable<3> constexpr omega {};  // angular velocity [rad/s]
    Variable<4> constexpr F {};      // horizontal force [N]
    
    auto constexpr p_dot = v;
    auto constexpr theta_dot = omega;
    auto const v_dot = (-l * m * sin(theta) * pow(omega, 2) + F + g * m * cos(theta) * sin(theta)) / (M + m - m * pow(cos(theta), 2));
    auto const omega_dot = (-l * m * cos(theta) * sin(theta) * pow(omega, 2) + F * cos(theta) + g * m * sin(theta) + M * g * sin(theta)) 
        / (l * (M + m - m * pow(cos(theta), 2)));

    auto arg = std::make_tuple(x[0], x[1], x[2], x[3], u[0]);

    dx[0] = eval(p_dot, arg);
    dx[1] = eval(theta_dot, arg);
    dx[2] = eval(v_dot, arg);
    dx[3] = eval(omega_dot, arg);

    ++vv;
}


double f1_dot_31(Real const x[4], Real const u[1])
{
    ++vv;

    auto const p = x[0];
    auto const theta = x[1];
    auto const v = x[2];
    auto const omega = x[3];
    auto const F = u[0];
    
    auto const a = -l * m * cos(theta) * sin(theta) * pow(omega, 2) + F * cos(theta) + g * m * sin(theta) + M * g * sin(theta);
    auto const b = l * (M + m - m * pow(cos(theta), 2));

    auto const a_dot = -l * m * (-sin(theta) * sin(theta) + cos(theta) * cos(theta)) * pow(omega, 2) 
        - F * sin(theta) + g * m * cos(theta) + M * g * cos(theta);
    auto const b_dot = l * m * 2 * cos(theta) * sin(theta);

    return (a_dot * b - a * b_dot) / (b * b);
}


double f2_dot_31(Real const x[4], Real const u[1])
{
    using tmpc::Variable;
    using tmpc::eval;

    ++vv;
    
    Variable<0> p {};      // horizontal displacement [m]
    Variable<1> theta {};  // angle with the vertical [rad]
    Variable<2> v {};      // horizontal velocity [m/s]
    Variable<3> omega {};  // angular velocity [rad/s]
    Variable<4> F {};      // horizontal force [N]
    
    auto omega_dot = (-l * m * cos(theta) * sin(theta) * pow(omega, 2) + F * cos(theta) + g * m * sin(theta) + M * g * sin(theta)) 
        / (l * (M + m - m * pow(cos(theta), 2)));

    auto domega_dot_dtheta = diff(omega_dot, theta);
    auto arg = std::make_tuple(x[0], x[1], x[2], x[3], u[0]);

    /*
    std::cout << "omega_dot = " << std::endl << omega_dot << std::endl;
    std::cout << "d(omega_dot)/d(theta) = " << std::endl << domega_dot_dtheta << std::endl;
    std::cout << "opCount = " << domega_dot_dtheta.operationCount() << std::endl;
    */

    return eval(domega_dot_dtheta, arg);
}


int main(int, char **)
{
    using clock = std::chrono::high_resolution_clock;
    using nanoseconds = std::chrono::duration<float, std::nano>;

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

    std::cout << "f1_dot_31(x, u): " << f1_dot_31(x, u) << std::endl;
    std::cout << "f2_dot_31(x, u): " << f2_dot_31(x, u) << std::endl;

    std::size_t N = 100000000;

    /*
    {
        std::cout << "Benchmarking f1()..." << std::endl;
        auto const f1_start = clock::now();
        for (std::size_t i = 0; i < N; ++i)
            f1(x, u, dx);
        auto const t1 = clock::now() - f1_start;
        std::cout << "f1() time per evaluation: " << nanoseconds(t1).count() / N << " ns" << std::endl;
    }

    {
        std::cout << "Benchmarking f2()..." << std::endl;
        auto const f2_start = clock::now();
        for (std::size_t i = 0; i < N; ++i)
            f2(x, u, dx);
        auto const t2 = clock::now() - f2_start;
        std::cout << "f2() time per evaluation: " << nanoseconds(t2).count() / N << " ns" << std::endl;
    }
    */

    {
        std::cout << "Benchmarking f1_dot_31()..." << std::endl;
        auto const f1_start = clock::now();
        for (std::size_t i = 0; i < N; ++i)
            result = f1_dot_31(x, u);
        auto const t1 = clock::now() - f1_start;
        std::cout << "f1_dot_31() time per evaluation: " << nanoseconds(t1).count() / N << " ns" << std::endl;
    }

    {
        std::cout << "Benchmarking f2_dot_31()..." << std::endl;
        auto const f2_start = clock::now();
        for (std::size_t i = 0; i < N; ++i)
            result = f2_dot_31(x, u);
        auto const t2 = clock::now() - f2_start;
        std::cout << "f2_dot_31() time per evaluation: " << nanoseconds(t2).count() / N << " ns" << std::endl;
    }

    {
        using tmpc::Variable;

        Variable<0> x;
        Variable<1> y;

        auto arg = std::make_tuple(3., 4.);

        Variable<0> p {};      // horizontal displacement [m]
        Variable<1> theta {};  // angle with the vertical [rad]
        Variable<2> v {};      // horizontal velocity [m/s]
        Variable<3> omega {};  // angular velocity [rad/s]
        Variable<4> F {};      // horizontal force [N]
        
        /*
        auto omega_dot = (-l * m * cos(theta) * sin(theta) * pow(omega, 2) + F * cos(theta) + g * m * sin(theta) + M * g * sin(theta)) 
            / (l * (M + m - m * pow(cos(theta), 2)));
        */
        auto omega_dot = log(theta);

        /*
        std::cout << "dx/dy = " << std::endl << diff(x, y) << std::endl;
        std::cout << "dx/dx = " << std::endl << diff(x, x) << std::endl;
        std::cout << "d(x+y)/dx = " << std::endl << diff(x + y, x) << std::endl;
        std::cout << "d(sin(x))/dx = " << std::endl << diff(sin(x), x) << std::endl;
        std::cout << "d(sin(x))/dy = " << std::endl << diff(sin(x), y) << std::endl;
        std::cout << "d(sin(x)+sin(y))/dy = " << std::endl << diff(sin(x) + sin(y), y) << std::endl;
        std::cout << "d(sin(x*y))/dx = " << std::endl << diff(sin(x*y), x) << std::endl;
        std::cout << "d(-x)/dx = " << std::endl << diff(-x, x) << std::endl;
        std::cout << "d(x/y)/dx = " << std::endl << diff(x / y, x) << std::endl;
        std::cout << "d(x/y)/dy = " << std::endl << diff(x / y, y) << std::endl;
        std::cout << "d(42*x)/dx = " << std::endl << diff(42 * x, x) << std::endl;
        std::cout << "d(pow(x,2))/dx = " << std::endl << diff(pow(x, 2), x) << std::endl;
        std::cout << "d(log(x))/dx = " << std::endl << diff(log(x), x) << std::endl;
        */

        /*
        {
            using tmpc::Pow;

            auto p = pow(x, 2);
            std::cout << "p = " << p << std::endl;

            auto a = Pow::diffL(p.left(), p.right()) * diff(p.left(), x);
            std::cout << "Pow::diffL(p.left(), p.right()) * diff(p.left(), x) = " 
                << a << std::endl;

            auto b = Pow::diffR(p.left(), p.right()) * diff(p.right(), x);
            std::cout << "Pow::diffR(p.left(), p.right()) * diff(p.right(), x) = "
                << b << std::endl;

            std::cout << "a + b = " << a + b << std::endl;

            std::cout << "diff(p,x) = " << diff(p, x) << std::endl;
        }
        */
    }
    
    return 0;
}