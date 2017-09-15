#pragma once

namespace tmpc
{
    class Rand<...>
    {
        void generate()
        {
            // Writing random data
            Rand<typename TestFixture::Matrix> rand_matrix;
            Rand<typename TestFixture::Vector> rand_vector;

            for (std::size_t i = 0; i < N; ++i)
            {
                auto const& sz = this->size_[i];
                auto const nx1 = i + 1 < N ? this->size_[i + 1].nx() : 0;
                auto& stage = this->qp_[i];

                stage.set_Q(Q[i] = rand_matrix.generate(sz.nx(), sz.nx()));
                stage.set_R(R[i] = rand_matrix.generate(sz.nu(), sz.nu()));
                stage.set_S(S[i] = rand_matrix.generate(sz.nx(), sz.nu()));
                stage.set_q(q[i] = rand_vector.generate(sz.nx()));
                stage.set_r(r[i] = rand_vector.generate(sz.nu()));

                stage.set_A(A[i] = rand_matrix.generate(nx1, sz.nx()));
                stage.set_B(B[i] = rand_matrix.generate(nx1, sz.nu()));
                stage.set_b(b[i] = rand_vector.generate(nx1));

                stage.set_lbx(x_min[i] = rand_vector.generate(sz.nx()));
                stage.set_ubx(x_max[i] = rand_vector.generate(sz.nx()));
                stage.set_lbu(u_min[i] = rand_vector.generate(sz.nu()));
                stage.set_ubu(u_max[i] = rand_vector.generate(sz.nu()));
                stage.set_lbd(d_min[i] = rand_vector.generate(sz.nc()));
                stage.set_ubd(d_max[i] = rand_vector.generate(sz.nc()));
            }
        }
    }
}