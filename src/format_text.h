// 3D World - colored/formatted text printing
// by Frank Gennari
// 05/14/2025
#pragma once

#include <string>

inline std::string format_red   (std::string const &s) {return "\033[31m" + s + "\033[0m";}
inline std::string format_green (std::string const &s) {return "\033[32m" + s + "\033[0m";}
inline std::string format_yellow(std::string const &s) {return "\033[33m" + s + "\033[0m";}
inline std::string format_blue  (std::string const &s) {return "\033[34m" + s + "\033[0m";}
inline std::string format_bold  (std::string const &s) {return "\033[1m"  + s + "\033[0m";}
inline std::string format_uline (std::string const &s) {return "\033[4m"  + s + "\033[0m";}
