#pragma once

#include <tmpc/core/detail/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/Traits.hpp>

#include <vector>


namespace tmpc
{
    template <typename Real, size_t NX, size_t NU, size_t NC = 0>
    class StaticOcpQp
    {
    public:
        StaticOcpQp(OcpGraph const& g)
        :   graph_(g)
        ,   vertexProperties_(num_vertices(g))
        ,   edgeProperties_(num_edges(g))
        {
        }


        auto size() const
        {
            return make_function_property_map<OcpVertexDescriptor>(
                [] (OcpVertexDescriptor v) 
                { 
                    return OcpSize {NX, NU, NC};
                }
            );
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto H()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::H_, vertexProperties());
        }


        auto H() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::H_, vertexProperties());
        }


        auto Q()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::Q_, vertexProperties());
        }


        auto Q() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::Q_, vertexProperties());
        }


        auto R()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto R() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto S()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto S() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto q()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto q() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto r()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto r() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto lx()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto lx() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto ux()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto ux() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto lu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto lu() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto uu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto uu() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto C()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto C() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto D()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto D() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto ld()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ld() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ud()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }


        auto ud() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }



		/*
        auto BA()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::BA_, edgeProperties());
        }


        auto BA() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::BA_, edgeProperties());
        }
		*/


        auto A()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto A() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto B()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto B() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto b()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


        auto b() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


    private:
        struct VertexPropertyBundle
        {
            // H = [ R,   S
            //       S^T, Q]
            blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NU + NX, NU + NX, blaze::columnMajor>> H_;

            decltype(blaze::submatrix<NU, NU, NX, NX>(H_)) Q_ = blaze::submatrix<NU, NU, NX, NX>(H_);
            decltype(blaze::submatrix<0, 0, NU, NU>(H_)) R_ = blaze::submatrix<0, 0, NU, NU>(H_);
            decltype(blaze::submatrix<0, NU, NU, NX>(H_)) S_ = blaze::submatrix<0, NU, NU, NX>(H_);
            blaze::StaticVector<Real, NX> q_;
            blaze::StaticVector<Real, NU> r_;

            blaze::StaticVector<Real, NX> lx_;
            blaze::StaticVector<Real, NX> ux_;
            blaze::StaticVector<Real, NU> lu_;
            blaze::StaticVector<Real, NU> uu_;

            blaze::StaticMatrix<Real, NC, NX, blaze::columnMajor> C_;
            blaze::StaticMatrix<Real, NC, NU, blaze::columnMajor> D_;
            blaze::StaticVector<Real, NC> ld_;
            blaze::StaticVector<Real, NC> ud_;
        };


        struct EdgePropertyBundle
        {
            blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor> A_;
            blaze::StaticMatrix<Real, NX, NU, blaze::columnMajor> B_;
			
			/*
            blaze::DynamicMatrix<Real, blaze::columnMajor> BA_;

            using Submatrix = decltype(blaze::submatrix(BA_, 0, 0, 1, 1));

            Submatrix A_;
            Submatrix B_;
			*/
			
            blaze::StaticVector<Real, NX> b_;
        };


        auto vertexProperties()
        {
            return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        }


        auto vertexProperties() const
        {
            return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        }


        auto edgeProperties()
        {
            return make_iterator_property_map(edgeProperties_.begin(), edgeIndex());
        }


        auto edgeProperties() const
        {
            return make_iterator_property_map(edgeProperties_.begin(), edgeIndex());
        }


        auto edgeIndex() const
        {
            return get(graph::edge_index, graph_);
        }


        OcpGraph graph_;
        std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<EdgePropertyBundle> edgeProperties_;
    };


    template <typename Real, size_t NX, size_t NU, size_t NC>
    struct RealOf<StaticOcpQp<Real, NX, NU, NC>>
    {
        using type = Real;
    };


    /// @brief Specialize copyQpProperties for StaticOcpQp<> source.
    ///
    /// This is necessary to avoid copying elements which are "logically" empty but physically non-empty.
    /// Such elements are S, R, r, lu, uu, D for leaf nodes.
    /// If we try to copy them to a dynamically-sized QP, we will get a runtime-error because of matrix/vector size mismatch.
    ///
    template <typename Real, size_t NX, size_t NU, size_t NC, typename QpDst>
	inline void copyQpProperties(StaticOcpQp<Real, NX, NU, NC> const& src, QpDst& dst)
	{
		auto const vert = graph::vertices(src.graph());
		copyProperty(src.Q(), dst.Q(), vert);
		copyProperty(src.q(), dst.q(), vert);
		copyProperty(src.lx(), dst.lx(), vert);
		copyProperty(src.ux(), dst.ux(), vert);
		copyProperty(src.C(), dst.C(), vert);
		copyProperty(src.ld(), dst.ld(), vert);
		copyProperty(src.ud(), dst.ud(), vert);

        // "A vertex of a path or tree is *internal* if it is not a leaf"
        // https://en.wikipedia.org/wiki/Glossary_of_graph_theory_terms
        //
        // Copy S, R, r, lu, uu, D only for internal vertices, since leaf vertices don't have u.
        for (auto v : vert)
        {
            if (out_degree(v, src.graph()) > 0)
            {
                put(dst.S(), v, get(src.S(), v));
                put(dst.R(), v, get(src.R(), v));
                put(dst.r(), v, get(src.r(), v));
                put(dst.lu(), v, get(src.lu(), v));
                put(dst.uu(), v, get(src.uu(), v));
                put(dst.D(), v, get(src.D(), v));
            }
        }        
		
		auto const edg = graph::edges(src.graph());
		copyProperty(src.A(), dst.A(), edg);
		copyProperty(src.B(), dst.B(), edg);
		copyProperty(src.b(), dst.b(), edg);
	}
}