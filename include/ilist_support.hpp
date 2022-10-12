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

    std::vector<std::size_t> operator()(dict_data_type c)
    {
        for (std::size_t r = 0; r < pfp.M.size(); r++)
        {
        
        }
    }
};

}

#endif /* end of include guard: _PFP_ILIST_SUPPORT_HH */
