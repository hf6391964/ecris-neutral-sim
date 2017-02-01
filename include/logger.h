#pragma once

#include <fstream>


class Logger {
    bool logging_ = false;
    std::ofstream logFile_;

    public:
        Logger(std::string filename);
        ~Logger();

        void setLogging(bool logging);

        template<typename T>
        friend Logger& operator<<(Logger &out, const T& value);
};

template<typename T>
Logger& operator<<(Logger &out, const T& value) {
    if (out.logging_ && out.logFile_.is_open()) {
        out.logFile_ << value;
    }

    return out;
}


extern Logger logger;

