#ifndef Logger_HPP
#define Logger_HPP 1

#pragma once
#include <string>
#include <fstream>
#include <mutex>

class Logger{
public:
    Logger();
    ~Logger();

    void setLogFile(const std::string& filename);
    void log(const std::string& message, bool forceprint=false);
    void setVerbose(bool enable);
    void setPrint(bool enable);

private:
    std::ofstream log_stream;
    bool log_enable = false;
    bool print_enable = false;
    std::mutex log_mutex;
};



#endif
