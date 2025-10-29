
#ifndef Parser_HPP
#define Parser_HPP 1

#include <unordered_map>
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <sstream>

using namespace std;
namespace fs = std::filesystem;

void parseArgs(int argc, char* argv[], unordered_map<string, string> &args);

#endif
