
#include "Logger.hpp"
#include <iostream>

Logger::Logger() {
    // 初期状態では標準出力（デバッグ用）
    log_stream = std::ofstream();  // 初期は閉じておく
}

Logger::~Logger() {
    if (log_stream.is_open()) {
        log_stream.close();
    }
}

// ログファイルを設定
void Logger::setLogFile(const std::string& filename) {
    if (log_stream.is_open()) {
        log_stream.close();  // 既に開いていた場合は閉じる
    }

    log_stream.open(filename, std::ios::out | std::ios::app);  // 追記モード
    if (!log_stream.is_open()) {
        std::cerr << "Failed to open log file: " << filename << std::endl;
    }
}

// ログを書き込む
void Logger::log(const std::string& message, bool forceprint) {
    std::lock_guard<std::mutex> guard(log_mutex);
    
	if (log_enable) {
		if (log_stream.is_open()) {
			log_stream << message;
            log_stream.flush();
		} else {
			std::cerr << "Log file is not open!" << std::endl;
		}
	}

	if(print_enable || forceprint){
		std::cout << message;
	}
}

void Logger::setVerbose(bool enable) {
    log_enable = enable;
}

void Logger::setPrint(bool enable) {
    print_enable = enable;
}












