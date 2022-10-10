//
//  unit_tests.cpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#include <iostream>

#include <utils.hpp>

#undef max
#undef min

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>


//------------------------------------------------------------------------------

struct listener : Catch::EventListenerBase
{
  using EventListenerBase::EventListenerBase;
  
  virtual void testCaseStarting(Catch::TestCaseInfo const& testInfo) override
  {
      std::cout << testInfo.tagsAsString() << " " << testInfo.name << std::endl;
  }
};
CATCH_REGISTER_LISTENER(listener)

//------------------------------------------------------------------------------

std::string testfiles_dir = "../tests/files";

//------------------------------------------------------------------------------

template <class char_type>
void read_fasta_file(const char *filename, std::vector<char_type>& v){
    FILE* fd;
    
    if ((fd = fopen(filename, "r")) == nullptr)
        spdlog::error("open() file " + std::string(filename) + " failed" );
    
    v.clear();

    char c;
    while (fread( &c, sizeof(char), 1,fd) == 1) {
        if(c == '>'){
            while(fread( &c, sizeof(char), 1,fd) == 1 && c != '\n');
        }else{
            v.push_back(c);
            while(fread( &c, sizeof(char), 1,fd) == 1 && c!= '\n') v.push_back(c);
        }
    }
    fclose(fd);
}

//------------------------------------------------------------------------------

#include <pfp.hpp>
#include <ra_support.hpp>
#include <sa_support.hpp>

//TEST_CASE( "pfp<uint8_t> RA to yeast", "PFP on yeast.fasta" )
//{
//    std::vector<uint8_t> yeast;
//    read_fasta_file(std::string(testfiles_dir + "/yeast.fasta").c_str(), yeast);
//
//    pfpds::pf_parsing<uint8_t> pfp(testfiles_dir + "/yeast.fasta", 10);
//    pfpds::pfp_ra_support<uint8_t> ra_support(pfp);
//
//    bool all_good = true;
//    for (std::size_t i = 0; i < yeast.size(); i++)
//    {
//        all_good = all_good and (yeast[i] == ra_support(i));
//        if (not all_good) { spdlog::error("mismatch at {} over {}", i, yeast.size()); break; }
//    }
//    REQUIRE(all_good);
//}
//
//TEST_CASE( "pfp<uint32_t> RA to yeast", "PFP on yeast.fasta.parse" )
//{
//    std::vector<uint32_t> yeast_parse;
//    pfpds::read_file(std::string(testfiles_dir + "/yeast.fasta.parse").c_str(), yeast_parse);
//
//    pfpds::pf_parsing<uint32_t> pfp(testfiles_dir + "/yeast.fasta.parse", 5);
//    pfpds::pfp_ra_support<uint32_t> ra_support(pfp);
//
//    bool all_good = true;
//    for (std::size_t i = 0; i < yeast_parse.size(); i++)
//    {
//        all_good = all_good and (yeast_parse[i] + 10 == ra_support(i));
//        if (not all_good) { spdlog::error("mismatch at {} over {}", i, yeast_parse.size()); break; }
//    }
//    REQUIRE(all_good);
//}
//
//#include <sdsl/suffix_arrays.hpp>
//TEST_CASE( "pfp<uint8_t> SA for yeast", "PFP on yeast.fasta" )
//{
//    std::size_t w = 10;
//
//    pfpds::pf_parsing<char> pfp(testfiles_dir + "/yeast.fasta", w);
//    pfpds::pfp_sa_support<char> sa_support(pfp);
//
//    // TEST sa_ds
//    std::vector<char> yeast;
//    read_fasta_file(std::string(testfiles_dir + "/yeast.fasta").c_str(), yeast);
//    yeast.insert(yeast.begin(), w - 1, 3);
//    yeast.insert(yeast.end(), w - 1, 5);
//    yeast.push_back(4);
//    yeast.push_back(0);
//
//    uint8_t num_bytes = 1;
//    sdsl::csa_wt<> csa;
//    sdsl::construct_im(csa, static_cast<const char *>(&yeast[0]), num_bytes);
//
//    bool all_good = true;
//    for (std::size_t i = 0; i < yeast.size(); ++i)
//    {
//        all_good = all_good and (sa_support(i) == ((csa[i] + (yeast.size()) - w + 1) % (yeast.size())));
//        if (not all_good) { break; }
//    }
//
//    REQUIRE(all_good);
//}
//
//TEST_CASE( "pfp<uint32_t> SA for yeast's parse", "PFP on yeast.fasta.parse" )
//{
//    std::size_t w = 5;
//
//    pfpds::pf_parsing<uint32_t> pfp(testfiles_dir + "/yeast.fasta.parse", w);
//    pfpds::pfp_sa_support<uint32_t> sa_support(pfp);
//
//    // TEST sa_ds
//    std::vector<uint32_t> yeast_parse;
//    pfpds::read_file(std::string(testfiles_dir + "/yeast.fasta.parse").c_str(), yeast_parse);
//    yeast_parse.push_back(0);
//
//    std::vector<uint_t> yeast_parse_SA;
//    yeast_parse_SA.resize(yeast_parse.size());
//    std::size_t alphabet_size = (*std::max_element(yeast_parse.begin(), yeast_parse.end())) + 1;
//    sacak_int(&yeast_parse[0],&yeast_parse_SA[0],yeast_parse.size(), alphabet_size);
//
//    bool all_good = true;
//    for (std::size_t i = 1; i < yeast_parse.size(); ++i)
//    {
//        all_good = all_good and (sa_support(i + w - 1) == yeast_parse_SA[i]);
//        if (not all_good) { break; }
//    }
//
//    REQUIRE(all_good);
//}


TEST_CASE( "pfp<uint8_t> from example", "PFP on example" )
{
    std::vector<char> text = {'A','C','G','T','T','C','G','C','A','A','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#',
                              'A','C','G','T','T','C','G','G','A','A','T','A','G','T','T','C','C','G','G','G','A','G','G','T','T','A','A','C','#',
                              'A','A','C','T','T','C','G','T','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','G','T','A','A','C','#',
                              'A','A','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','C','G','T','T','T','A','C','#',
                              'A','C','G','T','T','C','G','T','C','A','C','T','A','G','T','T','C','C','G','G','G','A','G','T','T','A','C','#',
                              'A','A','C','T','T','C','G','G','A','A','T','A','G','T','C','C','G','G','G','A','C','T','A','A','C','#',
                              'A','C','G','T','T','C','G','T','G','A','C','T','A','G','T','T','C','C','G','G','G','A','C','T','A','A','C','#',
                              'A','C','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#','#'
                              };
    std::vector<std::string> dict{
            "##AC", "AC##", "AC#AAC", "AC#AC", "ACCTAGT", "ACCTTCG", "ACG", "ACTAAC", "ACTAGT", "ACTTCG",
            "CGCAAC", "CGCAC", "CGGAATAGT", "CGGGAC", "CGGGAGGT", "CGGGAGT", "CGT", "GTAAC", "GTCAC", "GTCCG",
            "GTGAC", "GTTAAC", "GTTAC", "GTTCCG", "GTTCG", "GTTTAC",
    };
    std::vector<uint32_t> parse{0, 6, 16, 24, 10, 8, 19, 15, 22, 3, 6, 16, 24, 12, 23, 14, 21, 2, 9, 16, 18, 4, 19,
                                14, 17, 2, 9, 11, 4, 19, 13, 6, 16, 25, 3, 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9,
                                12, 19, 13, 7, 3, 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1, 0};


    std::vector<uint32_t> indices{0, 1, 2, 3, 4};
    std::vector<uint8_t> dict2 = {'#', '#', 'G', 'A', 'T', 'T', 'A', 'C', EndOfWord,
                                  'A', 'C', 'A', 'T', '#', EndOfWord,
                                  'A', 'G', 'A', 'T', 'A', '#', '#', EndOfWord,
                                  'T', '#', 'G', 'A', 'T', 'A', 'C', EndOfWord,
                                  'T', '#', 'G', 'A', 'T', 'T', 'A', 'G', EndOfWord, EndOfDict};
    sdsl::bit_vector b_d = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    const size_t n_phrases = 5;
    const std::vector<size_t> phrase_length{8, 5, 7, 7, 8};
    sdsl::bit_vector b_bwt = {1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,1};

    const std::vector<std::pair<size_t, std::pair<size_t, size_t>>> M{
            {2, {0, 0}}, {6, {2, 2}}, {7, {3, 3}}, {7, {4, 4}}, {3, {0, 0}},
            {2, {2, 3}}, {2, {4, 4}}, {3, {1, 1}}, {5, {0, 0}}, {4, {2, 2}},
            {5, {3, 3}}, {5, {4, 4}}, {4, {1, 1}}, {6, {0, 0}}, {5, {2, 2}},
            {6, {3, 3}}, {6, {4, 4}}, {2, {1, 1}}, {4, {0, 0}}, {3, {2, 3}},
            {3, {4, 4}}, {4, {3, 3}}, {4, {4, 4}}
    };
    const std::vector<uint32_t> bwt_p = {3, 1, 4, 5, 2, 2};

    // Extract the reverse of the phrases
    std::vector<std::pair<std::string,uint32_t>> rev_dict;
    int i = 0;
    int rank = 0;
    while(i < dict2.size()-1){
    std::string s;
    while(i < dict2.size()-1 && dict2[i] != EndOfWord) s.append(1,dict2[i++]);
    i++;
    reverse(s.begin(), s.end());
    rev_dict.push_back({s,rank++});
    }
    std::sort(rev_dict.begin(),rev_dict.end());

    std::vector<uint32_t> frequencies{0, 1, 2, 1, 1, 1};
    size_t w = 2;

    TEST_COUT << "Begin paper test" << std::endl;

    pf_parsing<> pf(dict2, parse, frequencies, w);
    pf_parsing<pfp_wt_custom> pf_c(dict2, parse, frequencies, w);
    TEST_COUT << "Pfp built" << std::endl;

    // TEST n
    EXPECT_EQ(pf.n, text.size());
    TEST_COUT << "Test n" << std::endl;

    // TEST n_phrases
    EXPECT_EQ(pf.dict.n_phrases(), n_phrases);
    TEST_COUT << "Test n_phrases" << std::endl;

    // TEST b_d
    for(size_t i = 0; i < pf.dict.b_d.size(); ++i){
    EXPECT_EQ(pf.dict.b_d[i], b_d[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_d" << std::endl;

    // TEST phrase_length
    for (size_t i = 0; i < phrase_length.size(); ++i)
    {
    EXPECT_EQ(pf.dict.length_of_phrase(i+1), phrase_length[i]);
    }
    TEST_COUT << "Test phrase_length" << std::endl;

    // TEST b_bwt
    for (size_t i = 0; i < pf.b_bwt.size(); ++i)
    {
    EXPECT_EQ(pf.b_bwt[i], b_bwt[i]) << "at position: " << i;
    }
    TEST_COUT << "Test b_bwt" << std::endl;

    // TEST M
    for (size_t i = 0; i < pf.M.size(); ++i)
    {
    EXPECT_EQ(pf.M[i].len, M[i].first) << "at position: " << i;
    EXPECT_EQ(pf.M[i].left, M[i].second.first) << "at position: " << i;
    EXPECT_EQ(pf.M[i].right, M[i].second.second) << "at position: " << i;
    }
    TEST_COUT << "Test M" << std::endl;

    // TEST colex_document_array_helper
    std::vector<uint_t> colex_id(pf.dict.n_phrases());
    std::vector<uint_t> inv_colex_id(pf.dict.n_phrases()); // I am using it as starting positions
    for (int i = 0, j = 0; i < pf.dict.d.size(); ++i)
    if (pf.dict.d[i + 1] == EndOfWord)
    {
    colex_id[j] = j;
    inv_colex_id[j++] = i;
    }

    pf.dict.colex_document_array_helper(inv_colex_id, colex_id, 0, pf.dict.n_phrases());
    for (size_t i = 0; i < colex_id.size(); ++i){
    EXPECT_EQ(colex_id[i], rev_dict[i].second) << "at position: " << i;
    }
    TEST_COUT << "Test colex_document_array_helper" << std::endl;

    std::vector<uint32_t> sa(text.size(), 0);
    std::iota(sa.begin(), sa.end(), 0);
    auto cyclic_sort = [&](const size_t a, const size_t b) {
        const auto max_cmp = text.size();
        for (size_t i = 0; i < max_cmp; ++i) {
            if (text[(a + i) % text.size()] != text[(b + i) % text.size()])
                return text[(a + i) % text.size()] < text[(b + i) % text.size()];
        }
    };
    std::sort(sa.begin(), sa.end(), cyclic_sort);

    pfp_lce_support<> lce_ds(pf);
    pfp_sa_support<> pf_sa(pf);

    pfp_lce_support<pfp_wt_custom> lce_ds_sdsl(pf_c);
    pfp_sa_support<pfp_wt_custom> pf_sa_sdsl(pf_c);

    // TEST BWT(P) - wavelet tree
    for (size_t i = 0; i < pf.w_wt.size(); ++i)
    {
    EXPECT_EQ(pf.w_wt[i], bwt_p[i]) << "at position: " << i;
    }
    TEST_COUT << "Test BWT(P)" << std::endl;

    // TEST BWT(P) - wavelet tree
    for (size_t i = 0; i < pf.w_wt.size(); ++i)
    {
    EXPECT_EQ(pf.w_wt[i], pf_c.w_wt[i]) << "at position: " << i;
    }
    TEST_COUT << "Test comparing WT custom vs. SDSL" << std::endl;

    for (size_t i = 0; i < pf_sa.size(); ++i) {
    EXPECT_EQ(pf_sa.sa(i), sa[i]) << "at position: " << i;
    }
    TEST_COUT << "Test PFP SA" << std::endl;

    for (size_t i = 0; i < pf_sa.size(); ++i) {
    EXPECT_EQ(pf_sa.sa(i), pf_sa_sdsl.sa(i)) << "at position: " << i;
    }
    TEST_COUT << "Test PFP comparing SA SDSL vs. custom" << std::endl;


    uint8_t num_bytes = 1;
    // build cst of the Text
    TEST_COUT << "Computing CSA of the text" << std::endl;
    std::vector<char> tmp_text = {'#', 'G', 'A', 'T', 'T', 'A', 'C', 'A', 'T', '#',
                                  'G', 'A', 'T', 'A', 'C', 'A', 'T', '#',
                                  'G', 'A', 'T', 'T', 'A', 'G', 'A', 'T', 'A', 0};
    sdsl::csa_wt<> csa;
    sdsl::construct_im(csa, static_cast<const char *>(&tmp_text[0]), num_bytes);
    for (size_t i = 0; i < pf_sa.size(); ++i) {
    EXPECT_EQ(pf_sa.sa(i), (csa[i] + (tmp_text.size()) - w + 1) % (tmp_text.size())) << "at position: " << i;
    }
    TEST_COUT << "Test PFP SA and cst" << std::endl;


    std::vector<uint32_t> isa(text.size(), 0);
    std::vector<uint32_t> lcp(text.size(), 0);
    for (size_t i = 0; i < sa.size(); ++i)
    isa[sa[i]] = i;

    LCP_array_cyclic_text(&text[0], isa, sa, text.size(), lcp);

    for (int i = 1; i < text.size() - 1; ++i)
    {
    EXPECT_EQ(lce_ds.lce(pf_sa.sa(i), pf_sa.sa(i + 1)), lcp[i + 1]) << "At positions: " << pf_sa.sa(i) << " " << pf_sa.sa(i + 1);
    }
    TEST_COUT << "Test LCE ds" << std::endl;
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    Catch::Session session;
    
    using namespace Catch::Clara;
    
    auto cli = session.cli() |
    Opt( testfiles_dir, "dir" ) ["--test-dir"] ("specify the directory containing the test dna sequences files");
    
    session.cli(cli);
    
    int returnCode = session.applyCommandLine(argc, argv);
    
    if( returnCode != 0 ) return returnCode;
    
    session.run();
}