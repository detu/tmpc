#pragma once


namespace tmpc
{
    ///
    /// \brief CRTP base for expression classes.
    ///
    template <typename Derived>
    class ExpressionBase
    {
    public:
        template <typename T>
        void evalTo(T& result) const
        {
            derived().implEvalTo(result);
        }


    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, an object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		ExpressionBase() = default;
		ExpressionBase(ExpressionBase const&) = default;


		Derived& derived()
		{
			return static_cast<Derived&>(*this);
		}


		Derived const& derived() const
		{
			return static_cast<Derived const&>(*this);	
        }
    };
}