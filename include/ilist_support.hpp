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

#include <utils.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <pfp.hpp>
#include <wt.hpp>

namespace pfpds
{

template<typename dict_data_type, class wt_t = pfp_wt_custom>
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
                    auto start = pfp.w_wt.range_select(pfp.M[r].left, pfp.M[r].right, 1);
                    std::size_t lex_rank = pfp.M[r].l_left;
                    
                    uint_t p_it = 0;
                    while (lex_rank <= pfp.M[r].l_right)
                    {
                        // colex subrange of the character we want
                        std::pair<std::size_t, std::size_t> colex_subrange = pfp.Q(r, c);
                        colex_subrange.second = colex_subrange.first + colex_subrange.second - 1; // Q stores left, length
                        
                        auto curr = pfp.w_wt.range_select(colex_subrange.first,colex_subrange.second, p_it + 1);
                        
                        uint_t tot_prev = 0;
                        for (dict_data_type col = 0; col < pfp.Q.elements_in_row(r); col++)
                        {
                            if (col != c)
                            {
                                dict_data_type active_col = pfp.Q.select_in_row(col);
                                std::pair<std::size_t, std::size_t> sr = pfp.Q(r, col);
                                sr.second = sr.first + sr.second - 1;
                                tot_prev += pfp.w_wt.range_count(sr.first, sr.second, curr);
                            }
                        }
                        
                        // print lex_rank
                        out_ilist.push_back(lex_rank + tot_prev);
                        
                        // update lex_rank
                        lex_rank += tot_prev + 1;
                    }
                    
                    
                    
                    while (partial_ilist.size() < num_of_c)
                    {
                        auto next = pfp.w_wt.range_select(pfp.M[r].left, pfp.M[r].right, 1);
                    }
                    
                    
                }
                // get the colex subrange.


                // first figure out the colex subrange for the character I want.

            }
        }
        
        return out_ilist;
    }
};

}

#endif /* end of include guard: _PFP_ILIST_SUPPORT_HH */
