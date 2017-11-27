#include <tmpc/print/Exception.hpp>

#include <exception>


namespace tmpc
{
    void print(std::ostream& os, const std::exception& e, int level)
    {
        os << std::string(level, ' ') << "exception: " << e.what() << '\n';
        try 
        {
            std::rethrow_if_nested(e);
        } 
        catch(const std::exception& e) 
        {
            print(os, e, level + 1);
        } 
        catch(...) 
        {
        }
    }
}