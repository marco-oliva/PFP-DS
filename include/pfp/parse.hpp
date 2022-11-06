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
    std::vector<uint_t> saP;
    std::vector<uint_t> isaP;
    std::vector<int_t> lcpP;
    sdsl::rmq_succinct_sct<> rmq_lcp_P;
    // sdsl::bit_vector b_p; // Starting position of each phrase in D
    // sdsl::bit_vector::rank_1_type rank_b_p;
    // sdsl::bit_vector::select_1_type select_b_p;
    bool saP_flag = false;
    bool isaP_flag = false;
    bool lcpP_flag = false;
    bool rmq_lcp_P_flag = false;
    
    size_t alphabet_size;
    
    typedef size_t size_type;
    
    // Default constructor for load
    parse() {}
    
    parse(  std::vector<uint32_t>& p_,
    size_t alphabet_size_,
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
    
    parse(  std::string filename,
    size_t alphabet_size_,
    bool saP_flag_ = true,
    bool isaP_flag_ = false,
    bool lcpP_flag_ = false,
    bool rmq_lcp_P_flag_ = false ):
    alphabet_size(alphabet_size_)
    {
        // Building dictionary from file
        std::string tmp_filename = filename + std::string(".parse");
        read_file(tmp_filename.c_str(), p);
        p.push_back(0); // this is the terminator for the sacak algorithm
        
        build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);
    }
    
    void build(bool saP_flag_, bool isaP_flag_, bool lcpP_flag_, bool rmq_lcp_P_flag_){
        
        if ((p.size() > (0x7FFFFFFF - 2)) and (sizeof(uint_t) == 4))
        {
            spdlog::error("Parse exceeds size allowed for 32 bits. Please use 64 bits executable.");
            std::exit(1);
        }
        else
        {
            if (sizeof(uint_t) == 4) { spdlog::info("Using 32 bits uint_t"); }
            else { spdlog::info("Using 64 bits uint_t"); }
        }
        
        
        // TODO: check if it has been already computed
        if(saP_flag_){
            saP.resize(p.size());
            // suffix array of the parsing.
            spdlog::info("Computing SA of the parsing");
            _elapsed_time(
            sacak_int(&p[0],&saP[0],p.size(),alphabet_size);
            );
        }
        
        assert(!isaP_flag_ || (saP_flag || saP_flag_) );
        if(isaP_flag_ && !isaP_flag){
            // inverse suffix array of the parsing.
            spdlog::info("Computing ISA of the parsing");
            _elapsed_time(
            {
            isaP.resize(p.size());
            for(int i = 0; i < saP.size(); ++i){
                isaP[saP[i]] = i;
            }
            }
            );
            
        }
        
        if(lcpP_flag_){
            lcpP.resize(p.size());
            // LCP array of the parsing.
            spdlog::info("Computing LCP of the parsing");
            _elapsed_time(
            LCP_array(&p[0], isaP, saP, p.size(), lcpP);
            );
        }
        
        
        assert(!rmq_lcp_P_flag_ || (lcpP_flag || lcpP_flag_));
        if(rmq_lcp_P_flag_ && ! rmq_lcp_P_flag){
            rmq_lcp_P_flag = true;
            spdlog::info("Computing RMQ over LCP of the parsing");
            // Compute the LCP rank of P
            _elapsed_time(
            rmq_lcp_P = sdsl::rmq_succinct_sct<>(&lcpP);
            );
        }
        
    }
    
    // Serialize to a stream.
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        
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
