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
#include <ilist_support.hpp>
#include <sdsl/suffix_arrays.hpp>

//TEST_CASE( "Q matix 1", "Sparse Matrix" )
//{
//    std::vector<std::pair<std::size_t, std::size_t>> table_test =
//    {
//        {0,0},{0,0},{1,0},{2,0},{3,0},{0,0},{0,0},
//        {4,0},{0,0},{1,0},{1,0},{3,0},{0,0},{0,0},
//        {0,0},{0,0},{0,0},{0,0},{3,0},{5,0},{0,0},
//        {0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
//        {0,0},{0,0},{1,0},{2,0},{3,0},{0,0},{1,0},
//    };
//    std::size_t rows = 5;
//    uint8_t columns = 7;
//
//
//    pfpds::pf_parsing<uint8_t>::Q_table table(rows, columns);
//    for (std::size_t i = 0; i < rows; i++)
//    {
//        for (uint8_t j = 0; j < columns; j++)
//        {
//            std::size_t pos = (i * columns) + j;
//            if (table_test[pos] != std::make_pair(0UL,0UL))
//            {
//                table.append(i, j, table_test[pos]);
//            }
//        }
//    }
//
//    bool all_good = true;
//    for (std::size_t i = 0; i < rows; i++)
//    {
//        for (uint8_t j = 0; j < columns; j++)
//        {
//            std::size_t pos = (i * columns) + j;
//            all_good = all_good and (table_test[pos] == table(i, j));
//        }
//    }
//    REQUIRE(all_good);
//
//    REQUIRE(table.elements_in_row(0) == 3);
//    REQUIRE(table.elements_in_row(1) == 4);
//    REQUIRE(table.elements_in_row(2) == 2);
//    REQUIRE(table.elements_in_row(3) == 0);
//    REQUIRE(table.elements_in_row(4) == 4);
//
//    REQUIRE(table.select_in_row(0, 0) == 2);
//    REQUIRE(table.select_in_row(0, 1) == 3);
//    REQUIRE(table.select_in_row(0, 2) == 4);
//
//    REQUIRE(table.select_in_row(1, 0) == 0);
//    REQUIRE(table.select_in_row(1, 1) == 2);
//    REQUIRE(table.select_in_row(1, 2) == 3);
//    REQUIRE(table.select_in_row(1, 3) == 4);
//
//
//}
//
//TEST_CASE( "Q matix 2", "Sparse Matrix" )
//{
//    std::vector<std::pair<std::size_t, std::size_t>> table_test =
//    {
//            {1,0},{0,0},{1,0},{2,0},{3,0},{0,0},{0,0},
//            {4,0},{1,0},{1,0},{1,0},{3,0},{1,0},{1,0},
//            {0,0},{0,0},{0,0},{0,0},{3,0},{0,0},{0,0},
//            {0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
//            {0,0},{0,0},{1,0},{2,0},{3,0},{0,0},{0,0},
//    };
//
//    std::size_t rows = 5;
//    uint8_t columns = 7;
//
//
//    pfpds::pf_parsing<uint8_t>::Q_table table(rows, columns);
//    for (std::size_t i = 0; i < rows; i++)
//    {
//        for (uint8_t j = 0; j < columns; j++)
//        {
//            std::size_t pos = (i * columns) + j;
//            if (table_test[pos] != std::make_pair(0UL,0UL))
//            {
//                table.append(i, j, table_test[pos]);
//            }
//        }
//    }
//
//    bool all_good = true;
//    for (std::size_t i = 0; i < rows; i++)
//    {
//        for (uint8_t j = 0; j < columns; j++)
//        {
//            std::size_t pos = (i * columns) + j;
//            all_good = all_good and (table_test[pos] == table(i, j));
//        }
//    }
//    REQUIRE(all_good);
//
//    REQUIRE(table.elements_in_row(0) == 4);
//    REQUIRE(table.elements_in_row(1) == 7);
//    REQUIRE(table.elements_in_row(2) == 1);
//    REQUIRE(table.elements_in_row(3) == 0);
//    REQUIRE(table.elements_in_row(4) == 3);
//
//    REQUIRE(table.select_in_row(2, 0) == 4);
//}
//
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
//TEST_CASE( "pfp<uint32_t> ilist for yeast's parse", "PFP on yeast.fasta.parse" )
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
    size_t w = 2;

    std::vector<char> text = {'A','C','G','T','T','C','G','C','A','A','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#',
                              'A','C','G','T','T','C','G','G','A','A','T','A','G','T','T','C','C','G','G','G','A','G','G','T','T','A','A','C','#',
                              'A','A','C','T','T','C','G','T','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','G','T','A','A','C','#',
                              'A','A','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','C','G','T','T','T','A','C','#',
                              'A','C','G','T','T','C','G','T','C','A','C','T','A','G','T','T','C','C','G','G','G','A','G','T','T','A','C','#',
                              'A','A','C','T','T','C','G','G','A','A','T','A','G','T','C','C','G','G','G','A','C','T','A','A','C','#',
                              'A','C','G','T','T','C','G','T','G','A','C','T','A','G','T','T','C','C','G','G','G','A','C','T','A','A','C','#',
                              'A','C','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#','#'
                              };

    std::vector<std::string> dict
    {
    "##AC", // 0
    "AC##", // 1
    "AC#AAC", // 2
    "AC#AC", // 3
    "ACCTAGT", // 4
    "ACCTTCG", // 5
    "ACG", // 6
    "ACTAAC", // 7
    "ACTAGT", // 8
    "ACTTCG", // 9
    "CGCAAC", // 10
    "CGCAC", // 11
    "CGGAATAGT", // 12
    "CGGGAC", // 13
    "CGGGAGGT", // 14
    "CGGGAGT", // 15
    "CGT", // 16
    "GTAAC", // 17
    "GTCAC", // 18
    "GTCCG", // 19
    "GTGAC", // 20
    "GTTAAC", // 21
    "GTTAC", // 22
    "GTTCCG", // 23
    "GTTCG", // 24
    "GTTTAC", // 25
    };

    std::vector<uint8_t> dict2;
    for (auto& phrase : dict)
    {
        for (auto& c : phrase)
        { dict2.push_back(c); }
        dict2.push_back(EndOfWord);
    }
    dict2.push_back(EndOfDict);

    std::vector<uint32_t> parse
    {
        0, 6, 16, 24, 10, 8, 19, 15, 22, 3, 6, 16, 24, 12, 23, 14, 21, 2, 9, 16, 18, 4, 19,
        14, 17, 2, 9, 11, 4, 19, 13, 6, 16, 25, 3, 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9,
        12, 19, 13, 7, 3, 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1
    };
    for (auto& p_id : parse) { p_id = p_id + 1; }

    std::vector<uint_t> frequencies(dict.size() + 1, 0);
    for (auto& p_id : parse) { frequencies[p_id] += 1; }

    std::vector<uint_t> colex_id = {1, 0, 3, 2, 10, 7, 17, 21, 11, 18, 13, 20, 22, 25, 6, 19, 23, 9, 5, 24, 15, 12, 8, 4, 16, 14};

    // build pfp
    parse.push_back(0);
    pfpds::pf_parsing<uint8_t> pfp(dict2, parse, frequencies, w);

    pfpds::pfp_ilist_support<uint8_t> ilist(pfp);

    REQUIRE(pfp.dict.colex_id == colex_id);

    std::string bwt = "CCCCC#CCCTTTTC###GGTATAATTACCC#####GGGAGCAAAGGTTTTTGGGTTTAAAAAAAAAAGGTGTTTTTTTTTAAA"
                      "ATTTTTCCCCCCCCTTTAAAAAAACCCAAAACAAGGGGGTGGGGGCCCCCGGGGGGGGCCCCCCCCAAGCCAAAAACGAAAAA"
                      "ACCCCCGTCCTTTTCACCCCACGGGTGGTGTGTTTTTTTTGGGGGTGGGGCCCGCGGG";

    std::vector<uint_t> ilist_dollar;
    for (std::size_t i = 0; i < bwt.size(); i++) { if (bwt[i] == '#') { ilist_dollar.push_back(i); } }
    REQUIRE(ilist_dollar == ilist('#'));

    std::vector<uint_t> ilist_A;
    for (std::size_t i = 0; i < bwt.size(); i++) { if (bwt[i] == 'A') { ilist_A.push_back(i); } }
    REQUIRE(ilist_A == ilist('A'));

    std::vector<uint_t> ilist_C;
    for (std::size_t i = 0; i < bwt.size(); i++) { if (bwt[i] == 'C') { ilist_C.push_back(i); } }
    REQUIRE(ilist_C == ilist('C'));

    std::vector<uint_t> ilist_G;
    for (std::size_t i = 0; i < bwt.size(); i++) { if (bwt[i] == 'G') { ilist_G.push_back(i); } }
    REQUIRE(ilist_G == ilist('G'));

    std::vector<uint_t> ilist_T;
    for (std::size_t i = 0; i < bwt.size(); i++) { if (bwt[i] == 'T') { ilist_T.push_back(i); } }
    REQUIRE(ilist_T == ilist('T'));
}

TEST_CASE( "pfp<uint32_t> SA for yeast's parse", "PFP on yeast.fasta.parse" )
{
    std::size_t w = 5;

    pfpds::pf_parsing<uint32_t> pfp(testfiles_dir + "/yeast.fasta.parse", w, 10);
    pfpds::pfp_ilist_support<uint32_t> ilist_support(pfp);

    // TEST sa_ds
    std::vector<uint32_t> yeast_parse;
    pfpds::read_file(std::string(testfiles_dir + "/yeast.fasta.parse").c_str(), yeast_parse);
    yeast_parse.push_back(0);

    std::vector<uint_t> yeast_parse_SA;
    yeast_parse_SA.resize(yeast_parse.size());
    std::size_t alphabet_size = (*std::max_element(yeast_parse.begin(), yeast_parse.end())) + 1;
    sacak_int(&yeast_parse[0],&yeast_parse_SA[0],yeast_parse.size(), alphabet_size);

    std::vector<uint32_t> yeast_parse_bwt(yeast_parse.size());
    for (std::size_t i = 0; i < yeast_parse.size(); i++)
    {
        if (yeast_parse_SA[i] == 0) { yeast_parse_bwt[i] = yeast_parse.back(); }
        else { yeast_parse_bwt[i] = yeast_parse[yeast_parse_SA[i] - 1]; }
    }

    std::map<uint32_t, std::vector<uint_t>> ilists;
    for(size_t i = 0; i < yeast_parse_bwt.size(); ++i)
    {
        ilists[yeast_parse_bwt[i]].push_back(i);
    }

    bool all_good = true;
    for (auto const& il : ilists)
    {
        uint32_t adjusted = il.first;
        if (adjusted > 0) { adjusted += pfp.shift; }
        std::vector<uint_t> pfp_ilist = ilist_support(adjusted);
        std::vector<uint_t> pfp_ilist_adjusted;
        for (auto e : pfp_ilist)
        {
            if (e > 1) { pfp_ilist_adjusted.push_back(e - (pfp.w - 1)); }
            else { pfp_ilist_adjusted.push_back(e); }
        }

        if (not pfp_ilist_adjusted.empty()) { all_good = all_good and (pfp_ilist_adjusted == il.second); }

        if (not all_good) { break; }
    }

    REQUIRE(all_good);
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