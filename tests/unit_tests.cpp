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

void read_fasta_file(const char *filename, std::vector<uint8_t>& v){
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

//#include <dictionary.hpp>
//#include <parse.hpp>
#include <pfp.hpp>
#include <ra_support.hpp>

TEST_CASE( "pfp<uint8_t> RA to yeast", "PFP on yeast.fasta" )
{
    std::vector<uint8_t> yeast;
    read_fasta_file(std::string(testfiles_dir + "/yeast.fasta").c_str(), yeast);

    pfpds::pf_parsing<uint8_t> pfp(testfiles_dir + "/yeast.fasta", 10);
    pfpds::pfp_ra_support<uint8_t> ra_support(pfp);

    bool all_good = true;
    for (std::size_t i = 0; i < yeast.size(); i++)
    {
        all_good = all_good and (yeast[i] == ra_support(i));
        if (not all_good) { spdlog::error("mismatch at {} over {}", i, yeast.size()); break; }
    }
    REQUIRE(all_good);
}

TEST_CASE( "pfp<uint32_t> RA to yeast", "PFP on yeast.fasta.parse" )
{
    std::vector<uint32_t> yeast_parse;
    pfpds::read_file(std::string(testfiles_dir + "/yeast.fasta.parse").c_str(), yeast_parse);
    
    pfpds::pf_parsing<uint32_t> pfp(testfiles_dir + "/yeast.fasta.parse", 5);
    pfpds::pfp_ra_support<uint32_t> ra_support(pfp);
    
    bool all_good = true;
    for (std::size_t i = 0; i < yeast_parse.size(); i++)
    {
        all_good = all_good and (yeast_parse[i] + 10 == ra_support(i));
        if (not all_good) { spdlog::error("mismatch at {} over {}", i, yeast_parse.size()); break; }
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