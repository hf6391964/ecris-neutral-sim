#include "logger.h"


Logger logger("sim_log.txt");


Logger::Logger(std::string filename) {
    logFile_.open(filename);
    if (!logFile_.is_open()) {
        throw std::invalid_argument("Could not open log file for writing");
    }
}

Logger::~Logger() {
    if (logFile_.is_open()) {
        logFile_.close();
    }
}

void Logger::setLogging(bool logging) {
    logging_ = logging;
}

