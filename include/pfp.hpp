/* pfp - prefix free parsing
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
   \file pfp.hpp
   \brief pfp.hpp define and build the prefix-free parsing data structures.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_HH
#define _PFP_HH

#include <utils.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <dictionary.hpp>
#include <parse.hpp>
#include <wt.hpp>

namespace pfpds
{

template<typename dict_data_type, class wt_t = pfp_wt_custom>
class pf_parsing{
public:
    struct M_entry_t{
        uint_t len;
        uint_t left; // left and right are the extremes of the range
        uint_t right;
        uint_t l_left;
        uint_t l_right;
    };

    class Q_table
    {
    private:
        std::vector<std::pair<std::size_t, std::size_t>> values;
        std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>> non_zero_tmp_queue;
        sdsl::sd_vector<> non_zero_positions;
        sdsl::rank_support_sd<> rank;
        sdsl::select_support_sd<> select;
        std::size_t rows;
        dict_data_type columns;
        
        sdsl::sd_vector_builder builder;
        bool built = false;

        std::size_t last_pos_set = 0;

    public:
        
        Q_table() {}
        
        Q_table(std::size_t rows, dict_data_type columns)
            : rows(rows), columns(columns)
        {}

        void init(std::size_t rows_, dict_data_type columns_)
        {
            this->rows = rows_;
            this->columns = columns_;
        }

        void append(std::size_t r, dict_data_type c, std::pair<std::size_t, std::size_t> v)
        {
            assert(not built);
            assert(not (v.first == 0 and v.second == 0));
            assert(r < rows and c < columns);

            std::size_t pos = (r * columns) + c;
            assert(last_pos_set == 0 or pos > last_pos_set);

            last_pos_set = pos;

            non_zero_tmp_queue.emplace_back(pos, v);
        }

        void build_static_structures()
        {
            if (built) { return; } built = true;

            this->builder = sdsl::sd_vector_builder(rows * columns, non_zero_tmp_queue.size());

            values.clear();
            for (auto& to_insert : non_zero_tmp_queue)
            {
                builder.set(to_insert.first);
                values.push_back(to_insert.second);
            }

            non_zero_positions = sdsl::sd_vector(builder);
            sdsl::util::init_support(rank, &non_zero_positions);
            sdsl::util::init_support(select, &non_zero_positions);
        }

        std::pair<std::size_t, std::size_t> operator()(std::size_t r, dict_data_type c)
        {
            if (not built) { build_static_structures(); }

            std::size_t pos = (r * columns) + c;
            assert(r < rows and c < columns);

            if (non_zero_positions[pos]) { return values[rank(pos)]; }
            else { return std::make_pair(0,0); }
        }
        
        dict_data_type elements_in_row(std::size_t r)
        {
            if (not built) { build_static_structures(); }
            
            std::size_t curr_row_pos = r * columns;
            std::size_t next_row_pos = (r + 1) * columns;
            return (rank(next_row_pos) - rank(curr_row_pos));
        }
        
        // 0 based
        dict_data_type select_in_row(std::size_t r, dict_data_type i)
        {
            if (not built) { build_static_structures(); }

            assert(i < elements_in_row(r));
            
            std::size_t curr_row_pos = r * columns;
            std::size_t prev = rank(curr_row_pos);
            
            return select(prev + i + 1) % columns;
        }

        bool non_zero(std::size_t r, dict_data_type c)
        {
            if (not built) { build_static_structures(); }

            std::size_t pos = (r * columns) + c;
            assert(r < rows and c < columns);
            return non_zero_positions[pos];
        }
        
        void print()
        {
            for (std::size_t i = 0; i < rows; i++)
            {
                std::size_t plotted = 0;
                for (std::size_t j = 0; j < columns; j++)
                {
                    
                    auto v = this->operator()(i,j);
                    if (v != std::make_pair(0UL,0UL)) { std::cout << char(j) << ": (" << v.first << "," << v.second << ")\t"; plotted += 1; }
                }
                if (plotted == 0) { std::cout << "empty line -----"; }
                std::cout << std::endl;
            }
        }
    };
    
    dictionary<dict_data_type> dict;
    parse pars;
    std::vector<uint_t> freq;
    size_t n; // Size of the text
    size_t w; // Size of the window
    
    sdsl::bit_vector b_bwt;
    sdsl::bit_vector::rank_1_type b_bwt_rank_1;
    sdsl::bit_vector::select_1_type b_bwt_select_1;
    std::vector<M_entry_t> M;
    
    Q_table Q;
    
    wt_t w_wt;
    
    sdsl::bit_vector b_p;
    sdsl::bit_vector::rank_1_type rank_b_p;
    sdsl::bit_vector::select_1_type select_b_p;
    
    typedef size_t size_type;
    
    // Default constructor for load
    pf_parsing() {}
    
    pf_parsing(std::vector<uint8_t> &d_,
    std::vector<uint32_t> &p_,
    std::vector<uint_t> &freq_,
    size_t w_) :
    dict(d_, w_),
    pars(p_, dict.n_phrases() + 1),
    freq(freq_),
    w(w_)
    {
        // Uploading the frequency file
        assert(freq[0] == 0);
        
        // Compute the length of the string;
        compute_n();
        
        spdlog::info("Computing b_p");
        _elapsed_time(compute_b_p());
        
        spdlog::info("Computing b_bwt and M of the parsing");
        _elapsed_time(build_b_bwt_and_M_and_Q());

        spdlog::info("Computing W of BWT(P)");
        _elapsed_time(build_W());

        // Clear unnecessary elements
        // clear_unnecessary_elements();
    }
    
    pf_parsing( std::string filename, size_t w_):
    dict(filename, w_),
    pars(filename,dict.n_phrases()+1),
    freq(dict.n_phrases() + 1,0),
    w(w_)
    {
        // Generating freq
        for (std::size_t i = 0; i < pars.p.size(); i++) { freq[pars.p[i]] += 1; }
        
        // Compute the length of the string;
        compute_n();
        
        // b_p(pfp.n,0);
        spdlog::info("Computing b_p");
        _elapsed_time(compute_b_p());
        
        spdlog::info("Computing b_bwt and M of the parsing");
        _elapsed_time(build_b_bwt_and_M_and_Q());

        spdlog::info("Computing W of BWT(P)");
        _elapsed_time(build_W());

        // Clear unnecessary elements
        clear_unnecessary_elements();
    }
    
    void compute_b_p() {
        // Build the bitvector storing the position of the beginning of each phrase.
        b_p.resize(this->n); // all should be initialized at false by sdsl
        for(size_t i = 0; i < b_p.size(); ++i)
            b_p[i] = false; // bug in resize
        b_p[0] = true; // phrase_0 becomes phrase 1
        
        size_t i = 0;
        
        for(int j = 0; j < pars.p.size()-2; ++j){ // -2 because the beginning of the last phrase is in position 0
            // p[i]: phrase_id
            assert(pars.p[j] != 0);
            // phrase_length: select_b_d(p[i]+1)-select_b_d(p[i]);
            i += dict.length_of_phrase(pars.p[j]) - w;
            b_p[i] = true;
        }
        
        // Build rank and select on Sp
        rank_b_p = sdsl::bit_vector::rank_1_type(&b_p);
        select_b_p = sdsl::bit_vector::select_1_type(&b_p);
    }
    
    void compute_n(){
        // Compute the length of the string;
        n = 0;
        for (int j = 0; j < pars.p.size() - 1; ++j)
        {
            // parse.p[j]: phrase_id
            assert(pars.p[j] != 0);
            n += dict.length_of_phrase(pars.p[j]) - w;
        }
        //n += w; // + w because n is the length including the last w markers
        //n += w - 1; // Changed after changind b_d in dict // -1 is for the first dollar + w because n is the length including the last w markers
    }

    void build_b_bwt_and_M_and_Q()
    {
        // Build the bitvector storing the position of the beginning of each phrase.
        b_bwt.resize(n);
        for (size_t i = 0; i < b_bwt.size(); ++i)
            b_bwt[i] = false; // bug in resize
        
        assert(dict.d[dict.saD[0]] == EndOfDict);
        size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t j = 0;
        
        size_t l_left  = 0;
        size_t l_right = 0;
        while (i < dict.saD.size())
        {
            size_t left = i;
            
            auto sn = dict.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = dict.daD[i] + 1;
            assert(phrase > 0 && phrase < freq.size()); // + 1 because daD is 0-based
            size_t suffix_length = dict.select_b_d(dict.rank_b_d(sn + 1) + 1) - sn - 1;
            if (dict.b_d[sn] || suffix_length < w)
            {
                ++i; // Skip
            }
            else
            {
                // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
                // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up
                b_bwt[j++] = true;
                j += freq[phrase] - 1; // the next bits are 0s
                i++;
                
                l_right += freq[phrase] - 1;
                
                if (i < dict.saD.size())
                {
                    auto new_sn = dict.saD[i];
                    auto new_phrase = dict.daD[i] + 1;
                    assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
                    size_t new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                    
                    while (i < dict.saD.size() && (dict.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
                    {
                        j += freq[new_phrase];
                        ++i;
                        
                        l_right += freq[new_phrase];
                        
                        if (i < dict.saD.size())
                        {
                            new_sn = dict.saD[i];
                            new_phrase = dict.daD[i] + 1;
                            assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
                            new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                        }
                    }
                }
                
                // Computing M
                size_t right = i - 1;
                M_entry_t m;
                m.len = suffix_length;
                m.left = dict.colex_daD[dict.rmq_colex_daD(left, right)];
                m.right = dict.colex_daD[dict.rMq_colex_daD(left, right)];
                m.l_left = l_left;
                m.l_right = l_right;
                
                M.push_back(m);
                
                l_left = l_right + 1;
                l_right = l_left;
            }
        }
        
        // rank & select support for b_bwt
        b_bwt_rank_1 = sdsl::bit_vector::rank_1_type(&b_bwt);
        b_bwt_select_1 = sdsl::bit_vector::select_1_type(&b_bwt);

        // now build Q
        Q.init(M.size(), dict.alphabet_size);

        for (std::size_t row = 0; row < M.size(); row++)
        {
            const M_entry_t& m = M[row];

            std::vector<std::pair<dict_data_type, std::size_t>> phrase_counts;
            for (std::size_t r = m.left; r <= m.right; r++)
            {
                auto phrase = dict.colex_id[r];
                std::size_t phrase_start = dict.select_b_d(phrase + 1);
                std::size_t phrase_length = dict.length_of_phrase(phrase + 1);
                dict_data_type c = dict.d[phrase_start + (phrase_length - m.len - 1)];

                if (phrase_counts.empty() or phrase_counts.back().first != c)
                {
                    phrase_counts.push_back(std::make_pair(c, 1));
                }
                else { phrase_counts.back().second += 1; }
            }
            
            // update matrix row
            std::size_t range_start = m.left;
            for (auto const& c_count : phrase_counts)
            {
                std::pair<std::size_t, std::size_t> q_element = std::make_pair(range_start, c_count.second);
                Q.append(row, c_count.first, q_element);
                range_start += c_count.second;
            }
        }
    }
    
    void build_W() {
        // create alphabet (phrases)
        std::vector<uint32_t> alphabet(dict.n_phrases());
        for (size_t i = 0; i < dict.n_phrases(); ++i) {
            alphabet[i] = dict.colex_id[i] + 1;
        }
        
        // create BWT(P)
        std::vector<uint32_t> bwt_p(pars.p.size() - 1, 0);
        for (size_t i = 1; i < pars.saP.size(); ++i) // TODO: shoud we count end symbol in this?
        {
            if (pars.saP[i] > 0)
                bwt_p[i - 1] = pars.p[pars.saP[i] - 1];
            else
                bwt_p[i - 1] = pars.p[pars.p.size() - 2]; // TODO: this should be -1 only if 0 stay in pars
        }
        
        w_wt.construct(alphabet, bwt_p);
    }
    
    void clear_unnecessary_elements(){
        dict.daD.clear();
        dict.colex_daD.clear();
        dict.colex_id.clear();
        dict.inv_colex_id.clear();
        // pars.saP.clear(); // It is needed in sa_support
        //    dict.rmq_colex_daD.clear();
        //    dict.rMq_colex_daD.clear();
    }
    
    // Serialize to a stream.
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        
        written_bytes += dict.serialize(out, child, "dictionary");
        written_bytes += pars.serialize(out, child, "parse");
        written_bytes += my_serialize(freq, out, child, "frequencies");
        written_bytes += sdsl::write_member(n, out, child, "n");
        written_bytes += sdsl::write_member(w, out, child, "w");
        written_bytes += b_bwt.serialize(out, child, "b_bwt");
        written_bytes += b_bwt_rank_1.serialize(out, child, "b_bwt_rank_1");
        written_bytes += b_bwt_select_1.serialize(out, child, "b_bwt_select_1");
        written_bytes += sdsl::serialize(M, out, child, "M");
        written_bytes += w_wt.serialize(out, child, "w_wt");
        written_bytes += b_p.serialize(out, child, "b_p");
        written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
        written_bytes += select_b_p.serialize(out, child, "select_b_p");
        // written_bytes += dict.serialize(out, child, "dictionary");
        // written_bytes += pars.serialize(out, child, "parse");
        // written_bytes += sdsl::serialize(freq, out, child, "frequencies");
        // written_bytes += sdsl::write_member(n, out, child, "n");
        // written_bytes += sdsl::write_member(w, out, child, "w");
        // written_bytes += b_bwt.serialize(out, child, "b_bwt");
        // written_bytes += b_bwt_rank_1.serialize(out, child, "b_bwt_rank_1");
        // written_bytes += b_bwt_select_1.serialize(out, child, "b_bwt_select_1");
        // written_bytes += sdsl::serialize(M, out, child, "M");
        // written_bytes += w_wt.serialize(out, child, "w_wt");
        // written_bytes += b_p.serialize(out, child, "b_p");
        // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
        // written_bytes += select_b_p.serialize(out, child, "select_b_p");
        
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }
    
    //! Load from a stream.
    void load(std::istream &in)
    {
        dict.load(in);
        pars.load(in);
        my_load(freq, in);
        sdsl::read_member(n, in);
        sdsl::read_member(w, in);
        b_bwt.load(in);
        b_bwt_rank_1.load(in, &b_bwt);
        b_bwt_select_1.load(in, &b_bwt);
        sdsl::load(M, in);
        w_wt.load(in);
        b_p.load(in);
        rank_b_p.load(in, &b_p);
        select_b_p.load(in, &b_p);
        // dict.load(in);
        // pars.load(in);
        // sdsl::load(freq, in);
        // sdsl::read_member(n, in);
        // sdsl::read_member(w, in);
        // b_bwt.load(in);
        // b_bwt_rank_1.load(in, &b_bwt);
        // b_bwt_select_1.load(in, &b_bwt);
        // sdsl::load(M, in);
        // w_wt.load(in);
        // b_p.load(in);
        // rank_b_p.load(in, &b_p);
        // select_b_p.load(in, &b_p);
    }
    
    std::string filesuffix() const
    {
        return ".pf.ds.other";
    }
    
    
};


// Specialization for pfp_wt_custom
template <>
std::string pf_parsing<pfp_wt_custom>::filesuffix() const
{
    return ".pf.ds";
}

// Specialization for pfp_wt_sdsl
template <>
std::string pf_parsing<pfp_wt_sdsl>::filesuffix() const
{
    return ".pf.wt_sdsl.ds";
}

// Specialization for pfp_wt_sdsl_2
template <>
std::string pf_parsing<pfp_wt_sdsl_2>::filesuffix() const
{
    return ".pf.wt_sdsl_2.ds";
}

using pf_parsing_custom = pf_parsing<pfp_wt_custom>;
using pf_parsing_sdsl = pf_parsing<pfp_wt_sdsl>;

}



#endif /* end of include guard: _PFP_HH */
