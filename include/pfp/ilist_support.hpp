/* pfp - prefix free parsing suffix array support
Copyright (C) 2020 Ondřej Cvacho, nope

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
\file ilist_support.hpp
\brief ilist_support.hpp define and build the prefix-free parsing suffix array support data structures.
\author Ondřej Cvacho, nope
\date 03/04/2020
*/


#ifndef _PFP_ILIST_SUPPORT_HH
#define _PFP_ILIST_SUPPORT_HH

#include "utils.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include "pfp.hpp"
#include "wt.hpp"

namespace pfpds
{

template<typename dict_data_type, class wt_t = pfp_wt_sdsl>
class pfp_ilist_support {
public:
    pf_parsing<dict_data_type, wt_t>& pfp;

    pfp_ilist_support(pf_parsing<dict_data_type, wt_t> & pfp_)
        : pfp(pfp_)
    { }

    size_t size() const { return pfp.n; }

    std::vector<uint_t> operator()(dict_data_type c)
    {
        std::vector<uint_t> out_ilist;

        if (c == 0) {out_ilist.push_back(1); return out_ilist; }

        for (std::size_t r = 0; r < pfp.M.size(); r++)
        {
            if (pfp.Q.non_zero(r, c))
            {
                // easy case
                if (pfp.Q.elements_in_row(r) == 1)
                {
                    for (uint_t e = pfp.M[r].l_left; e <= pfp.M[r].l_right; e++)
                    {
                        out_ilist.push_back(e);
                    }
                }
                // hard case, suffix is preceded by multiple characters
                else
                {
                    std::pair<std::size_t, std::size_t> colex_subrange = pfp.Q(r, c);
                    colex_subrange.second = colex_subrange.first + colex_subrange.second - 1;
                    std::vector<std::pair<std::size_t, std::size_t>> points = pfp.w_wt.range_search_2d(colex_subrange.first, colex_subrange.second, pfp.w_wt.size());
                    std::sort(points.begin(), points.end());
                    std::size_t num_of_c = points.size();
                    
                    std::size_t lex_rank = pfp.M[r].l_left;
                    
                    uint_t p_it = 0;
                    while (p_it < num_of_c)
                    {
                        std::pair<std::size_t, std::size_t> c_sr = pfp.Q(r, c);
                        c_sr.second = c_sr.first + c_sr.second - 1; // Q stores left, length
                        
                        auto curr = points[p_it].first;
                        
                        uint_t tot_prev = 0;
                        for (dict_data_type col = 0; col < pfp.Q.elements_in_row(r); col++)
                        {
                            dict_data_type active_col = pfp.Q.select_in_row(r,col);
                            if (active_col != c)
                            {
                                std::pair<std::size_t, std::size_t> sr = pfp.Q(r, active_col);
                                sr.second = sr.first + sr.second - 1;
                                tot_prev += pfp.w_wt.range_count(sr.first, sr.second, curr);
                            }
                        }
                        
                        // print lex_rank
                        out_ilist.push_back(lex_rank + p_it + tot_prev);
                        
                        p_it++;
                    }
                }
            }
        }
        
        return out_ilist;
    }
};

}

#endif /* end of include guard: _PFP_ILIST_SUPPORT_HH */
