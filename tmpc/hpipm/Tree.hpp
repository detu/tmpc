#pragma once

#include <hpipm_tree.h>

#include <tmpc/ocp/OcpTree.hpp>

#include <memory>


namespace tmpc :: hpipm
{
    struct Tree
    :   tree
    {
        Tree(OcpTree const& g);

    private:
        std::unique_ptr<node []> nodes_;
        std::unique_ptr<int []> kids_;
    };
}