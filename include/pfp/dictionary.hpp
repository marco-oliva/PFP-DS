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

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

namespace pfpds
{

template <typename data_type>
class dictionary
{

public:
    std::vector<data_type> d;
    std::vector<uint_t> saD;
    std::vector<uint_t> isaD;
    std::vector<int_da> daD;
    std::vector<int_t> lcpD;
    sdsl::rmq_succinct_sct<> rmq_lcp_D;
    sdsl::bit_vector b_d; // Starting position of each phrase in D
    sdsl::bit_vector::rank_1_type rank_b_d;
    sdsl::bit_vector::select_1_type select_b_d;
    std::vector<int_t> colex_daD;
    sdsl::rmq_succinct_sct<> rmq_colex_daD;
    sdsl::range_maximum_sct<>::type rMq_colex_daD;
    std::vector<uint_t> colex_id;
    std::vector<uint_t> inv_colex_id;
    std::size_t alphabet_size = 0;
    
    bool saD_flag = false;
    bool isaD_flag = false;
    bool daD_flag = false;
    bool lcpD_flag = false;
    bool rmq_lcp_D_flag = false;
    
    typedef size_t size_type;
    
    // default constructor for load.
    dictionary() {}
    
    dictionary( std::vector<data_type>& d_,
               size_t w,
               bool saD_flag_ = true,
               bool isaD_flag_ = true,
               bool daD_flag_ = true,
               bool lcpD_flag_ = true,
               bool rmq_lcp_D_flag_ = true ):
               d(d_)
    {
        build(saD_flag_, isaD_flag_, daD_flag_, lcpD_flag_, rmq_lcp_D_flag_);
        //assert(d[0] == Dollar);
    }
    
    dictionary(std::string filename,
               size_t w,
               bool saD_flag_ = true,
               bool isaD_flag_ = true,
               bool daD_flag_ = true,
               bool lcpD_flag_ = true,
               bool rmq_lcp_D_flag_ = true)
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
        
        build(saD_flag_, isaD_flag_, daD_flag_, lcpD_flag_, rmq_lcp_D_flag_);
    }
    
    inline size_t length_of_phrase(size_t id) {
        assert(id > 0);
        return select_b_d(id+1)-select_b_d(id) - 1; // to remove the EndOfWord
    }
    
    inline size_t n_phrases(){
        return rank_b_d(d.size()-1);
    }
    
    void build(bool saD_flag_, bool isaD_flag_, bool daD_flag_, bool lcpD_flag_, bool rmq_lcp_D_flag_){
        // Get alphabet size
        alphabet_size = (*std::max_element(d.begin(), d.end())) + 1;
        
        // Building the bitvector with a 1 in each starting position of each phrase in D
        b_d.resize(d.size());
        for(size_t i = 0; i < b_d.size(); ++i) b_d[i] = false; // bug in resize
        b_d[0] = true; // Mark the first phrase
        for(long long int i = 1; i < d.size(); ++i )
            b_d[i] = (d[i-1]==EndOfWord);
        b_d[d.size()-1] = true; // This is necessary to get the length of the last phrase
        
        rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
        select_b_d = sdsl::bit_vector::select_1_type(&b_d);
        
        assert(!rmq_lcp_D_flag_ || (lcpD_flag || lcpD_flag_));
        
        // TODO: check if it has been already computed
        if(saD_flag_ && daD_flag_ && lcpD_flag_){
            saD.resize(d.size());
            lcpD.resize(d.size());
            daD.resize(d.size());
            // suffix array, LCP array, and Document array of the dictionary.
            spdlog::info("Computing SA, LCP, and DA of dictionary");
            _elapsed_time(
            gsacak_templated<data_type>(&d[0],&saD[0],&lcpD[0],&daD[0],d.size(), alphabet_size)
            );
        }else if(saD_flag_ && lcpD_flag_){
            saD.resize(d.size());
            lcpD.resize(d.size());
            // suffix array and LCP array of the dictionary.
            spdlog::info("Computing SA, and LCP of dictionary");
            _elapsed_time(
            gsacak_templated<data_type>(&d[0],&saD[0],&lcpD[0],nullptr,d.size(), alphabet_size)
            );
        } else if(saD_flag_ && daD_flag_){
            saD.resize(d.size());
            daD.resize(d.size());
            // suffix array and LCP array of the dictionary.
            spdlog::info("Computing SA, and DA of dictionary");
            _elapsed_time(
            gsacak_templated<data_type>(&d[0],&saD[0],nullptr,&daD[0],d.size(), alphabet_size)
            );
        } else if(saD_flag_){
            saD.resize(d.size());
            // suffix array and LCP array of the dictionary.
            spdlog::info("Computing SA of dictionary");
            _elapsed_time(
            gsacak_templated<data_type>(&d[0],&saD[0],nullptr,nullptr,d.size(), alphabet_size)
            );
        }
        
        assert(!isaD_flag_ || (saD_flag || saD_flag_) );
        if(isaD_flag_ && !isaD_flag){
            // inverse suffix array of the dictionary.
            spdlog::info("Computing ISA of dictionary");
            _elapsed_time(
            {
                isaD.resize(d.size());
                for(int i = 0; i < saD.size(); ++i){
                    isaD[saD[i]] = i;
                }
            }
            );
        }
        
        assert(!rmq_lcp_D_flag_ || (lcpD_flag || lcpD_flag_));
        if(rmq_lcp_D_flag_ && ! rmq_lcp_D_flag){
            rmq_lcp_D_flag = true;
            
            spdlog::info("Computing RMQ over LCP of dictionary");
            // Compute the LCP rank of D
            _elapsed_time(
            rmq_lcp_D = sdsl::rmq_succinct_sct<>(&lcpD)
            );
        }
        
   
        // if(colex_daD_flag_){
        // co-lex document array of the dictionary.
        spdlog::info("Computing co-lex DA of dictionary");
        _elapsed_time(
        // {
        //   std::vector<uint_t>colex_id(n_phrases());
        //   std::vector<uint_t>inv_colex_id(n_phrases()); // I am using it as starting positions
        //   for(int i = 0, j = 0; i < d.size(); ++i )
        //     if(d[i+1]==EndOfWord){
        //       colex_id[j] = j;
        //       inv_colex_id[j++] = i;
        //     }
    
        //   colex_document_array_helper(inv_colex_id,colex_id,0,n_phrases());
    
        //   // computing inverse colex id
        //   for(int i = 0; i < colex_id.size(); ++i){
        //     inv_colex_id[colex_id[i]] = i;
        //   }
        //   colex_id.clear();
    
        //   colex_daD.resize(d.size());
        //   for(int i = 0; i < colex_daD.size(); ++i ){
        //     colex_daD[i]  = inv_colex_id[daD[i]];
        //   }
        // }
        {
            compute_colex_da();
            rmq_colex_daD = sdsl::rmq_succinct_sct<>(&colex_daD);
            rMq_colex_daD = sdsl::range_maximum_sct<>::type(&colex_daD);
        }
        );
        
        // }
    }
    
    void compute_colex_da()
    {
        colex_id.resize(n_phrases());
        inv_colex_id.resize(n_phrases());
        
//        std::vector<uint_t> inv_colex_id(n_phrases()); // I am using it as starting positions
//        for (long long int i = 0, j = 0; i < d.size(); ++i)
//            if (d[i + 1] == EndOfWord)
//            {
//                colex_id[j] = j;
//                inv_colex_id[j++] = i;
//            }
//
//        // buckets stores the begin and the end of each bucket.
//        std::queue<std::pair<long long int, long long int>> buckets;
//        // the first bucket is the whole array.
//        buckets.push({0,colex_id.size()});
//
//        // for each bucket
//        std::size_t iterations = 0;
//        while(!buckets.empty()){
//            auto bucket = buckets.front(); buckets.pop();
//            long long int start = bucket.first;
//            long long int end = bucket.second;
//            if ((start < end) && (end - start > 1))
//            {
//                std::vector<uint_t> count(alphabet_size, 0);
//                for (size_t i = start; i < end; ++i)
//                {
//                    count[d[inv_colex_id[i]]]++;
//                }
//
//                std::vector<uint_t> psum(alphabet_size, 0);
//                for (size_t i = 1; i < alphabet_size; ++i)
//                {
//                    psum[i] = psum[i - 1] + count[i - 1];
//                }
//
//                std::vector<uint_t> tmp(end - start, 0);
//                std::vector<uint_t> tmp_id(end - start, 0);
//                for (size_t i = start; i < end; ++i)
//                {
//                    size_t index = psum[d[inv_colex_id[i]]]++;
//                    tmp[index] = std::min(inv_colex_id[i] - 1, static_cast<uint_t>(d.size() - 1));
//                    tmp_id[index] = colex_id[i];
//                }
//
//                // Recursion
//                size_t tmp_start = 0;
//                for (size_t i = 0; i < alphabet_size; ++i)
//                {
//                    for (size_t j = 0; j < count[i]; ++j)
//                    {
//                        inv_colex_id[start + j] = tmp[tmp_start];
//                        colex_id[start + j] = tmp_id[tmp_start++];
//                    }
//                    end = start + count[i];
//                    if (i > EndOfWord){
//                        buckets.push({start, end});
//                    }
//                    start = end;
//                }
//            }
//
//        }
//
//        // computing inverse colex id
//        for (std::size_t i = 0; i < colex_id.size(); ++i)
//        {
//            inv_colex_id[colex_id[i]] = i;
//        }
//        colex_id.clear();
//
//        colex_daD.resize(d.size());
//        for (std::size_t i = 0; i < colex_daD.size(); ++i)
//        {
//            colex_daD[i] = inv_colex_id.at(daD[i] % inv_colex_id.size());
//        }

        // ---------- just sort the reversed phrases
        std::vector<std::pair<std::vector<data_type>,uint32_t>> rev_dict(n_phrases());
        std::size_t i = 0;
        std::size_t rank = 0;
        while(i < d.size()-1)
        {
            while((i < d.size()-1) and (d[i] != EndOfWord)) { rev_dict[rank].first.push_back(d[i++]); }
            i++;
            reverse(rev_dict[rank].first.begin(), rev_dict[rank].first.end());
            rev_dict[rank].second = rank;
            rank++;
        }
        std::sort(rev_dict.begin(),rev_dict.end());
        
        for (i = 0; i < colex_id.size(); i++) { colex_id[i] = rev_dict[i].second; }
        for (i = 0; i < colex_id.size(); i++) { inv_colex_id[colex_id[i]] = i; }
        
        // check
//        assert(inv_colex_id.size() == inv_colex_id_test.size());
//        for (i = 0; i < inv_colex_id_test.size(); i++)
//        {
//            if (inv_colex_id[i] != inv_colex_id_test[i])
//            {
//                std::cout << i << std::endl;
//            }
//            assert(inv_colex_id[i] == inv_colex_id_test[i]);
//        }
        
        colex_daD.resize(d.size());
        for (i = 0; i < colex_daD.size(); ++i)
        {
            colex_daD[i] = inv_colex_id.at(daD[i] % inv_colex_id.size());
        }
    }
    
    // Serialize to a stream.
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        
        written_bytes += my_serialize(d, out, child, "dictionary");
        written_bytes += my_serialize(saD, out, child, "saD");
        written_bytes += my_serialize(isaD, out, child, "isaD");
        written_bytes += my_serialize(daD, out, child, "daD");
        written_bytes += my_serialize(lcpD, out, child, "lcpD");
        written_bytes += rmq_lcp_D.serialize(out, child, "rmq_lcp_D");
        written_bytes += b_d.serialize(out, child, "b_d");
        written_bytes += rank_b_d.serialize(out, child, "rank_b_d");
        written_bytes += select_b_d.serialize(out, child, "select_b_d");
        written_bytes += my_serialize(colex_daD, out, child, "colex_daD");
        written_bytes += rmq_colex_daD.serialize(out, child, "rmq_colex_daD");
        written_bytes += rMq_colex_daD.serialize(out, child, "rMq_colex_daD");
        written_bytes += my_serialize(colex_id, out, child, "colex_id");
        written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");
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
        my_load(d, in);
        my_load(saD, in);
        my_load(isaD, in);
        my_load(daD, in);
        my_load(lcpD, in);
        rmq_lcp_D.load(in);
        b_d.load(in);
        rank_b_d.load(in, &b_d);
        select_b_d.load(in, &b_d);
        my_load(colex_daD, in);
        rmq_colex_daD.load(in);
        rMq_colex_daD.load(in);
        my_load(colex_id, in);
        sdsl::read_member(alphabet_size, in);
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
//void dictionary<char>::compute_colex_da(){
//    colex_id.resize(n_phrases());
//    std::vector<uint_t> inv_colex_id(n_phrases()); // I am using it as starting positions
//    for (int i = 0, j = 0; i < d.size(); ++i)
//        if (d[i + 1] == EndOfWord)
//        {
//            colex_id[j] = j;
//            inv_colex_id[j++] = i;
//        }
//
//    // buckets stores the begin and the end of each bucket.
//    std::queue<std::pair<int,int> > buckets;
//    // the first bucket is the whole array.
//    buckets.push({0,colex_id.size()});
//
//    // for each bucket
//    while(!buckets.empty()){
//        auto bucket = buckets.front(); buckets.pop();
//        int start = bucket.first;
//        int end = bucket.second;
//        if ((start < end) && (end - start > 1))
//        {
//            std::vector<uint32_t> count(256, 0);
//            for (size_t i = start; i < end; ++i)
//            {
//                count[d[inv_colex_id[i]]]++;
//            }
//
//            std::vector<uint32_t> psum(256, 0);
//            for (size_t i = 1; i < 256; ++i)
//            {
//                psum[i] = psum[i - 1] + count[i - 1];
//            }
//
//            std::vector<uint_t> tmp(end - start, 0);
//            std::vector<uint_t> tmp_id(end - start, 0);
//            for (size_t i = start; i < end; ++i)
//            {
//                size_t index = psum[d[inv_colex_id[i]]]++;
//                tmp[index] = std::min(inv_colex_id[i] - 1, static_cast<uint_t>(d.size() - 1));
//                tmp_id[index] = colex_id[i];
//            }
//
//            // Recursion
//            size_t tmp_start = 0;
//            for (size_t i = 0; i < 256; ++i)
//            {
//                for (size_t j = 0; j < count[i]; ++j)
//                {
//                    inv_colex_id[start + j] = tmp[tmp_start];
//                    colex_id[start + j] = tmp_id[tmp_start++];
//                }
//                end = start + count[i];
//                if (i > EndOfWord){
//                    buckets.push({start, end});
//                }
//                start = end;
//            }
//        }
//
//    }
//
//    // computing inverse colex id
//    for (int i = 0; i < colex_id.size(); ++i)
//    {
//        inv_colex_id[colex_id[i]] = i;
//    }
//    colex_id.clear();
//
//    colex_daD.resize(d.size());
//    for (int i = 0; i < colex_daD.size(); ++i)
//    {
//        colex_daD[i] = inv_colex_id[daD[i]];
//    }
//}

}

#endif /* end of include guard: _PFP_DICTIONARY_HH */
