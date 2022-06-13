/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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
