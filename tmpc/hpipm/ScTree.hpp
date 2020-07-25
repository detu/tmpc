#pragma once

#include <hpipm_scenario_tree.h>
#include <hpipm_tree.h>


namespace tmpc :: hpipm
{
    inline void create_sctree(int md, int Nr, int Nh, sctree& st, void * memory)
    {
        ::create_sctree(md, Nr, Nh, &st, memory);
    }
    
    
    inline void cast_sctree2tree(sctree const& st, tree& tt)
    {
        ::cast_sctree2tree(const_cast<sctree *>(&st), &tt);
    }
}