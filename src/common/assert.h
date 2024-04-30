#pragma once

#include <stdexcept>

#define ENSURE_OR_THROW(cond)                                                                      \
  do {                                                                                             \
    if (!(cond)) {                                                                                 \
      throw std::runtime_error(                                                                    \
              std::string(__FILE__) + ":" + std::to_string(__LINE__) + std::string(" ") + #cond); \
    }                                                                                              \
  } while (false)
