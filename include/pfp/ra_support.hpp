/*!
\file ra_support.hpp
\brief ra_support.hpp define and build the prefix-free parsing random access support data structures.
\author Ondřej Cvacho
\date 03/04/2020
*/


#ifndef _PFP_RA_SUPPORT_HH
#define _PFP_RA_SUPPORT_HH

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include "utils.hpp"
#include "pfp.hpp"
#include "wt.hpp"

namespace pfpds
{

template<typename dict_data_type, typename colex_comparator_type = std::less<dict_data_type>, class wt_t = pfp_wt_custom>
class pfp_ra_support {
public:
    const pf_parsing<dict_data_type, colex_comparator_type, wt_t>& pfp;
    
    pfp_ra_support(const pf_parsing<dict_data_type, colex_comparator_type, wt_t>& pfp_) : pfp(pfp_) {}
    
    long_type size() const { return pfp.n; }
    
    dict_data_type operator()(long_type i)
    {
        // get phrase rank
        long_type phrase_rank = pfp.rank_b_p.rank(i + pfp.w);

        assert(phrase_rank != 0);
        long_type phrase_offset = pfp.select_b_p(phrase_rank);
        uint32_t phrase_id = pfp.pars.p[phrase_rank - 1];
        
        long_type phrase_start = pfp.dict.select_b_d.select(phrase_id);
        
        return pfp.dict.d[phrase_start + (i + pfp.w - phrase_offset)];
    }
};

}



#endif /* end of include guard: _PFP_RA_SUPPORT_HH */
