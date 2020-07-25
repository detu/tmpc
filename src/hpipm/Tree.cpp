#include <tmpc/hpipm/Tree.hpp>


namespace tmpc :: hpipm
{
    Tree::Tree(OcpTree const& g)
    :   tree {}
    ,   nodes_ {new node[num_vertices(g)]}
    ,   kids_ {new int[num_vertices(g)]}
    {
        // The code of this function is based on create_sctree() from HPIPM:
        // https://github.com/giaf/hpipm/blob/master/tree_ocp_qp/scenario_tree.c
        
        Nn = num_vertices(g);
        memsize = 0;
        root = nodes_.get();
        kids = kids_.get();

        int tkids = 0;
        int idxkid = 0;

        // root
        node * node0 = root;
        node0->idx = 0;
        node0->stage = 0;
        node0->dad = -1;
        node0->real = -1;
        node0->idxkid = 0;

        // kids
        for (int idx = 0; idx < Nn; ++idx)
        {
            node0 = root + idx;
            auto const stage = node0->stage;
            auto const nkids = out_degree(idx, g);

            node0->nkids = nkids;
            if (nkids > 0)
            {
                node0->kids = kids + tkids;
                tkids += nkids;

                if (nkids > 1)
                {
                    for(int ii = 0; ii < nkids; ++ii)
                    {
                        ++idxkid;
                        node0->kids[ii] = idxkid;
                        node * node1 = root + idxkid;
                        node1->idx = idxkid;
                        node1->stage = stage + 1;
                        node1->dad = idx;
                        node1->real = ii;
                        node1->idxkid = ii;
                    }
                }
                else // nkids==1
                {
                    ++idxkid;
                    node0->kids[0] = idxkid;
                    node * node1 = root + idxkid;
                    node1->idx = idxkid;
                    node1->stage = stage + 1;
                    node1->dad = idx;
                    node1->real = node0->real;
                    node1->idxkid = 0;
                }
            }
        }
    }
}