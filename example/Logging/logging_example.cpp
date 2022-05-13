#include "bezierManipulation/src/utils/logger.hpp"


int main() {
  // Get the logger
  auto& logger = beziermanipulation::utils::Logger::Get();

  // Shorter version for the facilitated access
  using LoggerOutputOptions = beziermanipulation::utils::Logger::OutputLevel;
  logger.SetOutputLevel(15);
  logger.Warning("This is a warning");
  logger.Error("This is an error");
  logger.UserInfo("This is a user info");
  logger.UserInfo("Now enabling timestamps");
  logger.SetOutputLevel(
      {LoggerOutputOptions::all});
  logger.Warning("This is a warning");
  logger.Error("This is an error");
  logger.UserInfo("This is a user info");
  logger.Logging("This is a logging information");
  logger.ExtendedInformation("This is extended logging information");

  logger.TerminatingError<std::logic_error>("This is a terminating error! - Kills Program.");

  return 0;
}
