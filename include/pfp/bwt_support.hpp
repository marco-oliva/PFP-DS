//
// Created by marco on 10/27/22.
//

#ifndef PFP_DS_BWT_SUPPORT_HPP
#define PFP_DS_BWT_SUPPORT_HPP

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include "utils.hpp"
#include "pfp.hpp"
#include "sa_support.hpp"

namespace pfpds
{

    template<typename dict_data_type, class wt_t = pfp_wt_sdsl>
    class bwt_support {
    public:
        pf_parsing<dict_data_type, wt_t>& pfp;
        pfp_sa_support<dict_data_type, wt_t>& sa_support;


        bwt_support(pf_parsing<dict_data_type, wt_t> & pfp_)
                : pfp(pfp_), sa_support(pfp_)
        { }

        size_t size() const {
            return pfp.n;
        }

        dict_data_type operator()(std::size_t i)
        {
            auto sn = (sa_support(i) + pfp.w -1) % pfp.n; // suffix number
            auto p_i = pfp.rank_b_p(sn + 1);            // phrase number
            auto id_p_i = pfp.pars.p[p_i - 1];          // phrase_id of the phrase that i belongs to.
            size_t occ_in_p_i_in_D = pfp.dict.select_b_d(id_p_i) + (sn - pfp.select_b_p(p_i));
            return pfp.dict.d[occ_in_p_i_in_D];
        }
    };

}

#endif //PFP_DS_BWT_SUPPORT_HPP
