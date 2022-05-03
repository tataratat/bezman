#ifndef SRC_UTILS_LOGGER_HPP
#define SRC_UTILS_LOGGER_HPP

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace beziermanipulation::utils {

class Logger {
 public:
  // OutputLevel options
  enum class OutputLevel : unsigned int {
    nothing = 0,
    userinfo = 1,
    errors = 2,
    warnings = 4,
    logging = 8,
    logging_verbose = 16,
    time_stamp = 32,
    warning_file = 64,
    error_file = 128,
    log_file = 256
  };

  static Logger& GetLogger() {
    static Logger singleton_instance;
    singleton_instance.init();
    return singleton_instance;
  }

  void SetOutputLevel(unsigned int outputlevel) {
    outputlevel_ = outputlevel;
    init();
  }

  void SetOutputLevel(
      const std::initializer_list<OutputLevel>& output_options) {
    outputlevel_ = 0;
    for (auto i_option = output_options.begin();
         i_option != output_options.end(); ++i_option) {
      outputlevel_ |= static_cast<unsigned>(*i_option);
    }
    init();
  }

  void Warning(const std::string& warning_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::warnings) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_column)
                << ColorText("warning", Color::yellow) << "]" << GetTimeStamp()
                << " : " << ColorText(warning_text, Color::yellow) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::warning_file) & outputlevel_) {
      warning_file << "[" << std::setw(padding_first_column) << "warning"
                   << "]" << GetTimeStamp() << " : " << warning_text << "\n";
    }
#endif
  }

  void UserInfo(const std::string& info_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::userinfo) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_column)
                << ColorText("user info", Color::blue) << "]" << GetTimeStamp()
                << " : " << ColorText(info_text, Color::blue) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      log_file << "[" << std::setw(padding_first_column) << "user info"
               << "]" << GetTimeStamp() << " : " << info_text << "\n";
    }
#endif
  }

  void Logging(const std::string& log_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::logging) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_column)
                << ColorText("log", Color::shell_default) << "]"
                << GetTimeStamp() << " : "
                << ColorText(log_text, Color::shell_default) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      log_file << "[" << std::setw(padding_first_column) << "user info"
               << "]" << GetTimeStamp() << " : " << log_text << "\n";
    }
#endif
  }

  void ExtendedInformation(const std::string& log_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::logging_verbose) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_column)
                << ColorText("log", Color::green) << "]" << GetTimeStamp()
                << " : " << ColorText(log_text, Color::green) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      log_file << "[" << std::setw(padding_first_column) << "user info"
               << "]" << GetTimeStamp() << " : " << log_text << "\n";
    }
#endif
  }

  void Error(const std::string& error_text) {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::errors) & outputlevel_) {
      std::cout << "[" << std::setw(padding_first_column)
                << ColorText("error", Color::red) << "]" << GetTimeStamp()
                << " : " << ColorText(error_text, Color::red) << "\n";
    }
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      error_file << "[" << std::setw(padding_first_column) << "error"
                 << "]" << GetTimeStamp() << " : " << error_text << "\n";
    }
#endif
  }

  ~Logger() {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      if (is_init) {
        log_file.close();
      }
    }
    if (static_cast<unsigned>(OutputLevel::warning_file) & outputlevel_) {
      if (is_init) {
        warning_file.close();
      }
    }
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      if (is_init) {
        error_file.close();
      }
    }
#endif
  }

 private:
  /// Padding
  const unsigned int padding_first_column{18};

  /// File Streams
  std::ofstream log_file;
  std::ofstream warning_file;
  std::ofstream error_file;

  /// Initialized
  bool is_init{false};

  /// Initialize output streams
  void init() {
#ifdef ENABLE_LOGGING
    if (static_cast<unsigned>(OutputLevel::log_file) & outputlevel_) {
      if (!is_init) {
        log_file.open("log_file.log");
      }
    }
    if (static_cast<unsigned>(OutputLevel::warning_file) & outputlevel_) {
      if (!is_init) {
        warning_file.open("warning_file.log");
      }
    }
    if (static_cast<unsigned>(OutputLevel::error_file) & outputlevel_) {
      if (!is_init) {
        error_file.open("error_file.log");
      }
    }
#endif
  }

  /// Output Level
  unsigned int outputlevel_{39};

  /// Get the time stamp
  std::string GetTimeStamp() {
    if (static_cast<unsigned>(OutputLevel::time_stamp) & outputlevel_) {
      std::time_t result = std::time(nullptr);
      std::string time_stamp = (std::asctime(std::localtime(&result)));
      time_stamp.pop_back();
      return "[" + time_stamp + "]";
    } else {
      return std::string();
    }
  }

  // Available Colors
  enum Color { red, green, yellow, blue, white, cyan, magenta, shell_default };

  /// Color Text
  std::string ColorText(const std::string& text, Color color) {
    switch (color) {
      case Color::red:
        return std::string("\033[31m") + text + std::string("\033[0m");
      case Color::green:
        return std::string("\033[32m") + text + std::string("\033[0m");
      case Color::yellow:
        return std::string("\033[33m") + text + std::string("\033[0m");
      case Color::blue:
        return std::string("\033[34m") + text + std::string("\033[0m");
      case Color::magenta:
        return std::string("\033[35m") + text + std::string("\033[0m");
      case Color::cyan:
        return std::string("\033[36m") + text + std::string("\033[0m");
      case Color::white:
        return std::string("\033[37m") + text + std::string("\033[0m");
      case Color::shell_default:
        return std::string("\033[;0m") + text + std::string("\033[0m");
        return text;
      default:
        return text;
    }
  }

  Logger() {}
  Logger(Logger const&);          // Don't Implement.
  void operator=(Logger const&);  // Don't implement
};
}  // namespace beziermanipulation::utils
#endif  // SRC_UTILS_LOGGER_HPP
