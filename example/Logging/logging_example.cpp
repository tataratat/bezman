#include "bezman/src/utils/logger.hpp"

int main() {
  // Shorter version for the facilitated access
  using Logger = bezman::utils::Logger;
  using LoggerOutputOptions = bezman::utils::Logger::OutputLevel;
  Logger::SetOutputLevel(15);
  Logger::Warning("This is a warning");
  Logger::Error("This is an error");
  Logger::UserInfo("This is a user info");
  Logger::UserInfo("Now enabling timestamps");
  Logger::SetOutputLevel({LoggerOutputOptions::all});
  Logger::Warning("This is a warning");
  Logger::Error("This is an error");
  Logger::UserInfo("This is a user info");
  Logger::Logging("This is a logging information");
  Logger::ExtendedInformation("This is extended logging information");

  Logger::TerminatingError<std::logic_error>(
      "This is a terminating error! - Kills Program.");

  return 0;
}
