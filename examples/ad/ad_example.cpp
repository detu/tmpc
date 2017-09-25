#include <json.hpp>

#include <iostream>

using nlohmann::json;

template <typename Derived>
class ExpressionBase
{
public:
    Derived& derived()
    {
        return static_cast<Derived&>(*this);
    }

    Derived const& derived() const
    {
        return static_cast<Derived const&>(*this);
    }
};


template <std::size_t N_>
class Variable
:   public ExpressionBase<Variable<N_>>
{
public:
    static auto constexpr N = N_;

};


template <typename OP, typename A>
class UnaryOp
:   public ExpressionBase<UnaryOp<OP, A>>
{
public:
    using Arg = A;
    using Op = OP;

    UnaryOp(Arg const& arg)
    :   arg_(arg)
    {        
    }

    auto const& arg() const
    {
        return arg_;
    }

private:
    Arg arg_;
};


template <typename OP, typename L, typename R>
class BinaryOp
:   public ExpressionBase<BinaryOp<OP, L, R>>
{
public:
    using Left = L;
    using Right = R;
    using Op = OP;

    BinaryOp(Left const& left, Right const& right)
    :   left_(left)
    ,   right_(right)
    {        
    }

    auto const& left() const
    {
        return left_;
    }

    auto const& right() const
    {
        return right_;
    }

private:
    Left left_;
    Right right_;
};

struct Sin
{
    static std::string name()
    {
        return "sin";
    }
};

struct Cos
{
    static std::string name()
    {
        return "cos";
    }
};

struct Times
{
    static std::string name()
    {
        return "times";
    }
};

struct Div
{
    static std::string name()
    {
        return "Div";
    }
};

struct Plus
{
    static std::string name()
    {
        return "plus";
    }
};

struct Minus
{
    static std::string name()
    {
        return "minus";
    }
};

struct Pow
{
    static std::string name()
    {
        return "pow";
    }
};

template <typename Expr>
decltype(auto) sin(ExpressionBase<Expr> const& x)
{
    return UnaryOp<Sin, Expr>(x.derived());
}

template <typename Expr>
decltype(auto) cos(ExpressionBase<Expr> const& x)
{
    return UnaryOp<Cos, Expr>(x.derived());
}

template <typename T, typename Expr>
decltype(auto) operator*(T const& a, ExpressionBase<Expr> const& b)
{
    return BinaryOp<Times, T, Expr>(a, b.derived());
}

template <typename Expr, typename T>
decltype(auto) pow(ExpressionBase<Expr> const& a, T const& b)
{
    return BinaryOp<Pow, Expr, T>(a.derived(), b);
}

template <typename ExprA, typename ExprB>
decltype(auto) operator+(ExpressionBase<ExprA> const& a, ExpressionBase<ExprB> const& b)
{
    return BinaryOp<Plus, ExprA, ExprB>(a.derived(), b.derived());
}

template <typename T, typename Expr>
decltype(auto) operator-(T const& a, ExpressionBase<Expr> const& b)
{
    return BinaryOp<Minus, T, Expr>(a, b.derived());
}

template <typename ExprA, typename ExprB>
decltype(auto) operator/(ExpressionBase<ExprA> const& a, ExpressionBase<ExprB> const& b)
{
    return BinaryOp<Div, ExprA, ExprB>(a.derived(), b.derived());
}


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


int main(int, char **)
{
    using Real = double;
    
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
    
    Real const M = 1.;      // mass of the cart [kg]
    Real const m = 0.1;     // mass of the ball [kg]
    Real const g = 9.81;    // gravity constant [m/s^2]
    Real const l = 0.8;     // length of the rod [m]

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

    json j = omega_dot;
    //json j = {{"a", 10}};

    std::cout << j.dump(4) << std::endl;

	return 0;
}