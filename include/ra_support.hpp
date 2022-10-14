/*!
\file ra_support.hpp
\brief ra_support.hpp define and build the prefix-free parsing random access support data structures.
\author Ondřej Cvacho
\date 03/04/2020
*/


#ifndef _PFP_RA_SUPPORT_HH
#define _PFP_RA_SUPPORT_HH

#include <utils.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <pfp.hpp>
#include <wt.hpp>

namespace pfpds
{

template<typename dict_data_type, class wt_t = pfp_wt_sdsl>
class pfp_ra_support {
public:
    const pf_parsing<dict_data_type, wt_t>& pfp;
    
    pfp_ra_support(const pf_parsing<dict_data_type, wt_t>& pfp_) : pfp(pfp_) {}
    
    size_t size() const { return pfp.n; }
    
    dict_data_type operator()(std::size_t i)
    {
        // get phrase rank
        std::size_t phrase_rank = pfp.rank_b_p.rank(i + pfp.w);

        assert(phrase_rank != 0);
        std::size_t phrase_offset = pfp.select_b_p(phrase_rank);
        uint32_t phrase_id = pfp.pars.p[phrase_rank - 1];
        
        std::size_t phrase_start = pfp.dict.select_b_d.select(phrase_id);
        
        return pfp.dict.d[phrase_start + (i + pfp.w - phrase_offset)];
    }
};

}



#endif /* end of include guard: _PFP_RA_SUPPORT_HH */
