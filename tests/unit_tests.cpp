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

std::size_t p_global = 75;
std::size_t w_global = 20;

//------------------------------------------------------------------------------

#include <dictionary.hpp>

TEST_CASE( "dictionary<uint8_t> Build SA_D", "Dictionary" )
{
    pfpds::dictionary<uint8_t> dictionary(testfiles_dir + "/yeast.fasta", 10, true, false, false, false, false);
    REQUIRE(true);
}

TEST_CASE( "dictionary<uint32_t> Build SA_D", "Dictionary" )
{
    pfpds::dictionary<uint32_t> dictionary(testfiles_dir + "/yeast.fasta.parse", 5, true, false, false, false, false);
    REQUIRE(true);
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    Catch::Session session;
    
    using namespace Catch::Clara;
    
    auto cli = session.cli() |
    Opt( testfiles_dir, "dir" ) ["--test-dir"] ("specify the directory containing the test dna sequences files") |
    Opt( w_global, "int" ) ["-W"] ("specify w") |
    Opt( p_global, "int" ) ["-P"] ("specify p");
    
    session.cli(cli);
    
    int returnCode = session.applyCommandLine(argc, argv);
    
    if( returnCode != 0 ) return returnCode;
    
    spdlog::info("Tests running with w: {}\tp: {}", w_global, p_global);
    
    session.run();
}