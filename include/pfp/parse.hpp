/*!
   \file parse.hpp
   \brief parse.hpp define and build the prefix-free parse data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#ifndef _PFP_PARSE_HH
#define _PFP_PARSE_HH

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include "utils.hpp"

namespace pfpds
{

// TODO: Extend it to non-integer alphabets
class parse{
public:
    std::vector<uint32_t> p;
    std::vector<long_type> saP;
    std::vector<long_type> isaP;
    std::vector<long_type> lcpP;
    sdsl::rmq_succinct_sct<> rmq_lcp_P;
    // sdsl::bit_vector b_p; // Starting position of each phrase in D
    // sdsl::bit_vector::rank_1_type rank_b_p;
    // sdsl::bit_vector::select_1_type select_b_p;
    bool saP_flag = false;
    bool isaP_flag = false;
    bool lcpP_flag = false;
    bool rmq_lcp_P_flag = false;
    
    long_type alphabet_size;
    
    // Default constructor for load
    parse() {}
    
    parse(
    std::vector<uint32_t>& p_,
    long_type alphabet_size_,
    bool saP_flag_ = true,
    bool isaP_flag_ = false,
    bool lcpP_flag_ = false,
    bool rmq_lcp_P_flag_ = false ):
    p(p_),
    alphabet_size(alphabet_size_)
    {
        assert(p.back() == 0);
        build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);
    }
    
    parse(
    std::string filename,
    long_type alphabet_size_,
    bool saP_flag_ = true,
    bool isaP_flag_ = false,
    bool lcpP_flag_ = false,
    bool rmq_lcp_P_flag_ = false ):
    alphabet_size(alphabet_size_)
    {
        // Building parse from file
        std::string tmp_filename = filename + std::string(".parse");
        read_file(tmp_filename.c_str(), p);
        p.push_back(0); // this is the terminator for the sacak algorithm
        
        build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);
    }
    
    void build(bool saP_flag_, bool isaP_flag_, bool lcpP_flag_, bool rmq_lcp_P_flag_)
    {
        // SA
        if(saP_flag_)
        {
            saP.resize(p.size());
            spdlog::info("Using 8 bytes for SA");
            sacak_int(&p[0], &saP[0], p.size(), alphabet_size);
            saP_flag = true;
        }
        
        // ISA
        if(isaP_flag_)
        {
            assert(saP_flag);
            spdlog::info("Computing ISA of the parsing");
            isaP.resize(p.size());
            for (long_type i = 0; i < saP.size(); i++) { isaP[saP[i]] = i; }
            isaP_flag = true;
        }
    
        // LCP
        if(lcpP_flag_)
        {
            spdlog::info("Computing LCP of the parsing");
            lcpP.resize(p.size());
            LCP_array(&p[0], isaP, saP, p.size(), lcpP);
        }
        
        // RMQ over LCP
        assert(!rmq_lcp_P_flag_ || (lcpP_flag || lcpP_flag_));
        if(rmq_lcp_P_flag_)
        {
            spdlog::info("Computing RMQ over LCP of the parsing");
            rmq_lcp_P = sdsl::rmq_succinct_sct<>(&lcpP);
            rmq_lcp_P_flag = true;
        }
        
    }
    
    // Serialize to a stream.
    long_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        long_type written_bytes = 0;
        
        written_bytes += my_serialize(p, out, child, "parse");
        written_bytes += my_serialize(saP, out, child, "saP");
        written_bytes += my_serialize(isaP, out, child, "isaP");
        written_bytes += my_serialize(lcpP, out, child, "lcpP");
        written_bytes += rmq_lcp_P.serialize(out, child, "rmq_lcp_P");
        // written_bytes += b_p.serialize(out, child, "b_p");
        // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
        // written_bytes += select_b_p.serialize(out, child, "select_b_p");
        written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");
        // written_bytes += sdsl::serialize(p, out, child, "parse");
        // written_bytes += sdsl::serialize(saP, out, child, "saP");
        // written_bytes += sdsl::serialize(isaP, out, child, "isaP");
        // written_bytes += sdsl::serialize(lcpP, out, child, "lcpP");
        // written_bytes += rmq_lcp_P.serialize(out, child, "rmq_lcp_P");
        // // written_bytes += b_p.serialize(out, child, "b_p");
        // // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
        // // written_bytes += select_b_p.serialize(out, child, "select_b_p");
        // written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");
        
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }
    
    //! Load from a stream.
    void load(std::istream &in)
    {
        my_load(p, in);
        my_load(saP, in);
        my_load(isaP, in);
        my_load(lcpP, in);
        rmq_lcp_P.load(in);
        // b_p.load(in);
        // rank_b_p.load(in);
        // select_b_p.load(in);
        sdsl::read_member(alphabet_size, in);
        // sdsl::load(p, in);
        // sdsl::load(saP, in);
        // sdsl::load(isaP, in);
        // sdsl::load(lcpP, in);
        // rmq_lcp_P.load(in);
        // // b_p.load(in);
        // // rank_b_p.load(in);
        // // select_b_p.load(in);
        // sdsl::read_member(alphabet_size, in);
    }
    
};

}

#endif /* end of include guard: _PFP_PARSE_HH */
