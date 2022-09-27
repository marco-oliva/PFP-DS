//
//  utils.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <string>
#include <string_view>
#include <vector>
#include <mutex>
#include <set>
#include <limits>
#include <functional>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <list>
#include <forward_list>
#include <unordered_map>
#include <cstring>

#include <unistd.h>
#include <math.h>

#include <spdlog/spdlog.h>

namespace pfpds
{

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter

typedef std::uint64_t long_type;
typedef std::uint32_t short_type;
typedef std::int64_t long_signed_type;
typedef std::int32_t short_signed_type;
typedef uint8_t  char_type;

#ifdef PFP_LONG_TYPE
//pragma message("Using 64 bit parse")
typedef long_type  size_type;
typedef long_signed_type signed_type;
#else
//#pragma message("Using 32 bit parse")
typedef short_type size_type;
typedef short_signed_type signed_type;
#endif
typedef long_type  hash_type;

constexpr std::size_t KILOBYTE      = 1024;
constexpr std::size_t MEGABYTE      = KILOBYTE * KILOBYTE;
constexpr std::size_t GIGABYTE      = KILOBYTE * MEGABYTE;

constexpr double KILOBYTE_DOUBLE    = 1024.0;
constexpr double MILLION_DOUBLE     = 1000000.0;
constexpr double MEGABYTE_DOUBLE    = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
constexpr double GIGABYTE_DOUBLE    = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

constexpr std::size_t MILLION       = 1000000;
constexpr std::size_t BILLION       = 1000 * MILLION;

namespace EXT
{

const std::string PARSE = ".parse";
const std::string DICT = ".dict";
const std::string DICT_COMPRESSED = ".dicz";
const std::string DICT_COMPRESSED_LENGTHS = ".dicz.len";
const std::string OCC = ".occ";
const std::string LAST = ".last";
const std::string SAI = ".sai";
const std::string N_PARSE = ".aup.parse";
const std::string N_DICT = ".aup.dict";
const std::string N_OCC = ".aup.occ";
const std::string N_LAST = ".aup.last";
const std::string N_SAI = ".aup.sai";
const std::string N_DICT_COMPRESSED = ".aup.dicz";
const std::string N_DICT_COMPRESSED_LENGTHS = ".aup.dicz.len";

}


//------------------------------------------------------------------------------

inline double
inMegabytes(std::size_t bytes)
{
    return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(std::size_t bytes)
{
    return bytes / GIGABYTE_DOUBLE;
}

inline std::size_t pid()
{
#ifdef MSVC_COMPILER
    return _getpid();
#else
    return getpid();
#endif
}

//------------------------------------------------------------------------------

/*
  Temporary file names have the pattern "prefix_hostname_pid_counter", where
  - prefix is given as an argument to getName();
  - hostname is the name of the host;
  - pid is the process id; and
  - counter is a running counter starting from 0.
  The generated names are stored until the file is deleted with remove(). All
  remaining temporary files are deleted when the program exits (normally or
  with std::exit()).
  TempFile is thread-safe.
*/

namespace TempFile
{
extern const std::string DEFAULT_TEMP_DIR;
extern std::string temp_dir;

void setDirectory(const std::string& directory);
std::string getName(const std::string& name_part);
void remove(std::string& filename);  // Also clears the filename.
}

void truncate_file(std::string file_name, std::size_t new_size_in_bytes);

bool is_gzipped(std::ifstream& in);

//------------------------------------------------------------------------------

// Disk space statistics
namespace DiskWrites
{
void update(std::size_t num_of_bytes);
};

//------------------------------------------------------------------------------

} // end namespace vcfbwt

#endif //utils_hpp