#pragma once

#include "common/task.h"

namespace ocm {

namespace heuristic {

  namespace local {

    Positions Optimize(const Task& task, const Positions& positions);

  }

}// namespace heuristic

}// namespace ocm