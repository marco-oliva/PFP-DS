/* pfp-dictionary - prefix free parsing dictionary
    Copyright (C) 2020 Massimiliano Rossi

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
   \file dictionary.hpp
   \brief dictionary.hpp define and build the prefix-free dictionary data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_DICTIONARY_HH
#define _PFP_DICTIONARY_HH

#include <queue>

#include "utils.hpp"
#undef max

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

namespace pfpds
{

template <typename data_type, class colex_comparator_type = std::less<data_type>>
class dictionary
{

public:
    std::vector<data_type> d;
    std::vector<long_type> saD;
    sdsl::int_vector<0> isaD;
    sdsl::int_vector<0> daD;
    sdsl::int_vector<0> lcpD;
    sdsl::rmq_succinct_sct<> rmq_lcp_D;
    sdsl::bit_vector b_d; // Starting position of each phrase in D
    sdsl::bit_vector::rank_1_type rank_b_d;
    sdsl::bit_vector::select_1_type select_b_d;
    sdsl::int_vector<0> colex_daD;
    sdsl::rmq_succinct_sct<> rmq_colex_daD;
    sdsl::range_maximum_sct<>::type rMq_colex_daD;
    sdsl::int_vector<0> colex_id;
    sdsl::int_vector<0> inv_colex_id;
    long_type alphabet_size = 0;
    
    bool saD_flag = false;
    bool isaD_flag = false;
    bool daD_flag = false;
    bool lcpD_flag = false;
    bool rmq_lcp_D_flag = false;
    bool colex_id_flag = false;
    bool colex_daD_flag = false;
    
    long_type w;
    
    colex_comparator_type& colex_comparator;
    
    // default constructor for load.
    dictionary() {}
    
    dictionary(std::vector<data_type>& d_,
               long_type w,
               colex_comparator_type& colex_comparator,
               bool saD_flag_ = true,
               bool isaD_flag_ = true,
               bool daD_flag_ = true,
               bool lcpD_flag_ = true,
               bool rmq_lcp_D_flag_ = true,
               bool colex_id_flag_ = true,
               bool colex_daD_flag = true):
               d(d_), w(w), colex_comparator(colex_comparator)
    {
        build(saD_flag_, isaD_flag_, daD_flag_, lcpD_flag_, rmq_lcp_D_flag_, colex_id_flag_, colex_daD_flag);
        //assert(d[0] == Dollar);
    }
    
    dictionary(std::string filename,
               long_type w,
               colex_comparator_type& colex_comparator,
               bool saD_flag_ = true,
               bool isaD_flag_ = true,
               bool daD_flag_ = true,
               bool lcpD_flag_ = true,
               bool rmq_lcp_D_flag_ = true,
               bool colex_id_flag_ = true,
               bool colex_daD_flag = true):
              w(w), colex_comparator(colex_comparator)
    {
        // Building dictionary from file
        std::string tmp_filename = filename + std::string(".dict");
        read_file(tmp_filename.c_str(), d);
        assert(d[0] == Dollar);
        // Prepending w dollars to d
        // 1. Count how many dollars there are
        int i = 0;
        int n_dollars = 0;
        while(i < d.size() && d[i++] == Dollar)
            ++n_dollars;
        std::vector<data_type> dollars(w-n_dollars,Dollar);
        d.insert(d.begin(), dollars.begin(),dollars.end());
        
        build(saD_flag_, isaD_flag_, daD_flag_, lcpD_flag_, rmq_lcp_D_flag_, colex_id_flag_, colex_daD_flag);
    }
    
    inline long_type length_of_phrase(long_type id) {
        assert(id > 0);
        return select_b_d(id+1)-select_b_d(id) - 1; // to remove the EndOfWord
    }
    
    inline long_type n_phrases(){
        return rank_b_d(d.size()-1);
    }
    
    void build(bool saD_flag_, bool isaD_flag_, bool daD_flag_, bool lcpD_flag_, bool rmq_lcp_D_flag_, bool colex_id_flag_, bool colex_daD_flag_){
      
        // Get alphabet size
        alphabet_size = (*std::max_element(d.begin(), d.end())) + 1;
        
        // Building the bitvector with a 1 in each starting position of each phrase in D
        // Also compute length of longhest phrase
        b_d.resize(d.size());
        long_type max_phrase_length = 0;
        for(long_type i = 0; i < b_d.size(); ++i) { b_d[i] = false; } // bug in resize
        b_d[0] = true; // Mark the first phrase
        long_type prev_start = 0;
        for(long_type i = 1; i < d.size(); ++i )
        {
            if (d[i-1] == EndOfWord)
            {
                b_d[i] = true;
                max_phrase_length = std::max(max_phrase_length, i - prev_start - 2);
                prev_start = i;
            }
        }
        b_d[d.size()-1] = true; // This is necessary to get the length of the last phrase
        
        rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
        select_b_d = sdsl::bit_vector::select_1_type(&b_d);
        
        // SA
        if(saD_flag_)
        {
            _elapsed_time(
            spdlog::info("Using 8 bytes for SA of the dictionary");
            saD.resize(d.size());
            gsacak_templated<data_type>(&d[0], &saD[0], d.size(), alphabet_size);
            saD_flag = true;
            );
        }
    
        // DA
        if(daD_flag_)
        {
            _elapsed_time(
            assert(saD_flag);
            long_type bytes_daD = 0;
            long_type ps = n_phrases(); assert(ps != 0);
            while (ps != 0) { ps >>= 8; bytes_daD++; }
            spdlog::info("Using {} bytes for DA of the dictionary", bytes_daD);
            daD = sdsl::int_vector<>(d.size(), 0, bytes_daD * 8);
            for (long_type i = 0; i < saD.size(); i++)
            {
                long_type out = rank_b_d.rank(saD[i]);
                if (b_d[saD[i]]) { out += 1; }
                daD[i] = out - 1;
            }
            daD_flag = true;
            );
        }
        
        // ISA
        if (isaD_flag_)
        {
            _elapsed_time(
            assert(saD_flag);
            long_type bytes_isaD = 0;
            long_type max_sa = d.size() + 1;
            while (max_sa != 0) { max_sa >>= 8; bytes_isaD++; }
            spdlog::info("Using {} bytes for ISA of the dictionary", bytes_isaD);
            isaD = sdsl::int_vector<>(d.size(), 0, bytes_isaD * 8);
            for (long_type i = 0; i < saD.size(); i++) { isaD[saD[i]] = i; }
            isaD_flag = true;
            );
        }
        
        // LCP
        if (lcpD_flag_)
        {
            _elapsed_time(
            assert(saD_flag and isaD_flag);
            long_type bytes_lcpD = 0;
            long_type mpl = max_phrase_length;
            while (mpl != 0) { mpl >>= 8; bytes_lcpD++; }
            spdlog::info("Using {} bytes for LCP of the dictionary", bytes_lcpD);
            lcpD = sdsl::int_vector<>(d.size(), 0, bytes_lcpD * 8);
    
            // Kasai et al. LCP construction algorithm
            lcpD[0]  = 0;
            long_type l = 0;
            for (long_type i = 0; i < lcpD.size(); ++i)
            {
                // if i is the last character LCP is not defined
                long_type k = isaD[i];
                if(k > 0)
                {
                    long_type j = saD[k-1];
                    // I find the longest common prefix of the i-th suffix and the j-th suffix.
                    while(d[i+l] == d[j+l] and d[i + l] != EndOfWord) { l++; }
                    // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
                    lcpD[k] = l;
                    if(l>0) { l--; }
                }
            }
            lcpD_flag = true;
            );
        }
        
        // RMQ over LCP
        if(rmq_lcp_D_flag_)
        {
            _elapsed_time(
            assert(lcpD_flag);
            spdlog::info("Computing RMQ over LCP of dictionary");
            rmq_lcp_D = sdsl::rmq_succinct_sct<>(&lcpD);
            rmq_lcp_D_flag = true;
            );
        }
    
        // co-lex document array of the dictionary.
        if(colex_daD_flag_ or colex_id_flag_)
        {
            _elapsed_time(
            assert(daD_flag);
            
            // allocating space for colex DA
            if (colex_daD_flag_)
            {
                long_type bytes_colex_daD = 0;
                long_type ps = n_phrases(); assert(ps != 0);
                while (ps != 0) { ps >>= 8; bytes_colex_daD++; }
                spdlog::info("Using {} bytes for colex DA of the dictionary", bytes_colex_daD);
                colex_daD = sdsl::int_vector<>(d.size(), 0, bytes_colex_daD * 8);
            }
            
            // allocating space for colex_id and inv_colex_id
            if (colex_daD_flag_ or colex_id_flag_)
            {
                long_type bytes_inv_colex_id = 0;
                long_type ps = n_phrases(); assert(ps != 0);
                while (ps != 0) { ps >>= 8; bytes_inv_colex_id++; }
                spdlog::info("Using {} bytes for colex id and inverse colex id array of the dictionary", bytes_inv_colex_id);
                colex_id = sdsl::int_vector<>(n_phrases(), 0, bytes_inv_colex_id * 8);
                inv_colex_id = sdsl::int_vector<>(n_phrases(), 0, bytes_inv_colex_id * 8);
            }
            
            compute_colex_da(colex_id_flag_, colex_daD_flag_);
            
            if (colex_daD_flag_)
            {
                rmq_colex_daD = sdsl::rmq_succinct_sct<>(&colex_daD);
                rMq_colex_daD = sdsl::range_maximum_sct<>::type(&colex_daD);
            }
            
            if (colex_daD_flag_) { colex_daD_flag = true; }
            if (colex_id_flag_) { colex_id_flag = true; }
            );
        }
    }
    
    void compute_colex_da(bool colex_id_flag_, bool colex_daD_flag_)
    {

        // ---------- just sort the reversed phrases
        std::vector<std::pair<std::vector<data_type>,uint32_t>> rev_dict(n_phrases());
        long_type i = 0;
        long_type rank = 0;
        while(i < d.size()-1)
        {
            while((i < d.size()-1) and (d[i] != EndOfWord)) { rev_dict[rank].first.emplace_back(d[i++]); }
            i++;
            reverse(rev_dict[rank].first.begin(), rev_dict[rank].first.end());
            rev_dict[rank].second = rank;
            rank++;
        }
        
        std::vector<std::pair<std::vector<data_type>,uint32_t>> rev_dict_2(rev_dict);
        
        std::sort(rev_dict.begin(),rev_dict.end(),
                  [this] (std::pair<std::vector<data_type>,uint32_t>& l, std::pair<std::vector<data_type>,uint32_t>& r)
                  {
                      std::pair<typename std::vector<data_type>::iterator, typename std::vector<data_type>::iterator> mismatch;
                      if (l.first.size() < r.first.size())
                      {
                          mismatch = std::mismatch(l.first.begin(), l.first.end(), r.first.begin());
                          if (mismatch.first == l.first.end())
                          {
                              return true;
                          }
                          return this->colex_comparator(*(mismatch.first), *(mismatch.second));
                      }
                      else
                      {
                          mismatch = std::mismatch(r.first.begin(), r.first.end(), l.first.begin());
                          if (mismatch.first == r.first.end())
                          {
                              return false;
                          }
                          return this->colex_comparator(*(mismatch.second), *(mismatch.first));
                      }
                  } );
        
        for (i = 0; i < colex_id.size(); i++) { colex_id[i] = rev_dict[i].second; }
        for (i = 0; i < colex_id.size(); i++) { inv_colex_id[colex_id[i]] = i; }
        
        for (i = 0; i < colex_daD.size(); ++i)
        {
            colex_daD[i] = inv_colex_id[daD[i] % inv_colex_id.size()];
        }
    }
    
    // Serialize to a stream.
    long_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        long_type written_bytes = 0;
        
//        written_bytes += my_serialize(d, out, child, "dictionary");
//        written_bytes += my_serialize(saD, out, child, "saD");
//        written_bytes += my_serialize(isaD, out, child, "isaD");
//        written_bytes += my_serialize(daD, out, child, "daD");
//        written_bytes += my_serialize(lcpD, out, child, "lcpD");
//        written_bytes += rmq_lcp_D.serialize(out, child, "rmq_lcp_D");
//        written_bytes += b_d.serialize(out, child, "b_d");
//        written_bytes += rank_b_d.serialize(out, child, "rank_b_d");
//        written_bytes += select_b_d.serialize(out, child, "select_b_d");
//        written_bytes += my_serialize(colex_daD, out, child, "colex_daD");
//        written_bytes += rmq_colex_daD.serialize(out, child, "rmq_colex_daD");
//        written_bytes += rMq_colex_daD.serialize(out, child, "rMq_colex_daD");
//        written_bytes += my_serialize(colex_id, out, child, "colex_id");
//        written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");
        // written_bytes += sdsl::serialize(d, out, child, "dictionary");
        // written_bytes += sdsl::serialize(saD, out, child, "saD");
        // written_bytes += sdsl::serialize(isaD, out, child, "isaD");
        // written_bytes += sdsl::serialize(daD, out, child, "daD");
        // written_bytes += sdsl::serialize(lcpD, out, child, "lcpD");
        // written_bytes += rmq_lcp_D.serialize(out, child, "rmq_lcp_D");
        // written_bytes += b_d.serialize(out, child, "b_d");
        // written_bytes += rank_b_d.serialize(out, child, "rank_b_d");
        // written_bytes += select_b_d.serialize(out, child, "select_b_d");
        // written_bytes += sdsl::serialize(colex_daD, out, child, "colex_daD");
        // written_bytes += rmq_colex_daD.serialize(out, child, "rmq_colex_daD");
        // written_bytes += rMq_colex_daD.serialize(out, child, "rMq_colex_daD");
        // written_bytes += sdsl::serialize(colex_id, out, child, "colex_id");
        
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
        
    }
    
    //! Load from a stream.
    void load(std::istream &in)
    {
//        my_load(d, in);
//        my_load(saD, in);
//        my_load(isaD, in);
//        my_load(daD, in);
//        my_load(lcpD, in);
//        rmq_lcp_D.load(in);
//        b_d.load(in);
//        rank_b_d.load(in, &b_d);
//        select_b_d.load(in, &b_d);
//        my_load(colex_daD, in);
//        rmq_colex_daD.load(in);
//        rMq_colex_daD.load(in);
//        my_load(colex_id, in);
//        sdsl::read_member(alphabet_size, in);
        // sdsl::load(d, in);
        // sdsl::load(saD, in);
        // sdsl::load(isaD, in);
        // sdsl::load(daD, in);
        // sdsl::load(lcpD, in);
        // rmq_lcp_D.load(in);
        // b_d.load(in);
        // rank_b_d.load(in, &b_d);
        // select_b_d.load(in, &b_d);
        // sdsl::load(colex_daD, in);
        // rmq_colex_daD.load(in);
        // rMq_colex_daD.load(in);
        // sdsl::load(colex_id, in);
    }
};


//template <>
//void dictionary<uint8_t>::compute_colex_da(bool colex_id_flag_, bool colex_daD_flag_){
//
//    for (long_type i = 0, j = 0; i < d.size(); ++i)
//    {
//        if (d[i + 1] == EndOfWord)
//        {
//            colex_id[j] = j;
//            inv_colex_id[j++] = i;
//        }
//    }
//
//    // buckets stores the start and the end of each bucket.
//    std::queue<std::pair<long_type,long_type>> buckets;
//    // the first bucket is the whole array.
//    buckets.push({0,colex_id.size()});
//
//    // for each bucket
//    while(not buckets.empty())
//    {
//        auto bucket = buckets.front(); buckets.pop();
//        long_type start = bucket.first;
//        long_type  end = bucket.second;
//        if ((start < end) && (end - start > 1))
//        {
//            std::vector<uint32_t> count(256, 0);
//            for (long_type i = start; i < end; ++i)
//            {
//                count[d[inv_colex_id[i]]]++;
//            }
//
//            std::vector<uint32_t> psum(256, 0);
//            for (long_type i = 1; i < 256; ++i)
//            {
//                psum[i] = psum[i - 1] + count[i - 1];
//            }
//
//            std::vector<long_type > tmp(end - start, 0);
//            std::vector<long_type > tmp_id(end - start, 0);
//            for (long_type i = start; i < end; ++i)
//            {
//                long_type index = psum[d[inv_colex_id[i]]]++;
//                tmp[index] = std::min(inv_colex_id[i] - 1, (long_type) d.size() - 1);
//                tmp_id[index] = colex_id[i];
//            }
//
//            // Recursion
//            long_type tmp_start = 0;
//            for (long_type  i = 0; i < 256; ++i)
//            {
//                for (long_type  j = 0; j < count[i]; ++j)
//                {
//                    inv_colex_id[start + j] = tmp[tmp_start];
//                    colex_id[start + j] = tmp_id[tmp_start++];
//                }
//                end = start + count[i];
//                if (i > EndOfWord) { buckets.push({start, end}); }
//                start = end;
//            }
//        }
//
//    }
//
//    // computing inverse colex id
//    for (long_type  i = 0; i < colex_id.size(); ++i) { inv_colex_id[colex_id[i]] = i; }
//
//    if (colex_daD_flag_)
//    {
//        for (long_type  i = 0; i < colex_daD.size(); ++i)
//        {
//            colex_daD[i] = inv_colex_id[daD[i] % inv_colex_id.size()];
//        }
//    }
//
//}

}

#endif /* end of include guard: _PFP_DICTIONARY_HH */
