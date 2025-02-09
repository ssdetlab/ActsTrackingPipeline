#pragma once

#include <optional>

#include <tbb/parallel_for.h>
#include <tbb/queuing_mutex.h>
#include <tbb/task_arena.h>

namespace tbbWrap {

/// enableTBB keeps a record of whether we are multi-threaded (nthreads!=1) or
/// not. This is set once in task_arena and stored globally.
/// This means that enableTBB(nthreads) itself is not thread-safe. That should
/// be fine because the task_arena is initialised before spawning any threads.
/// If multi-threading is ever enabled, then it is not disabled.
static bool enableTBB(int nthreads = -99) {
  static bool setting = false;
  if (nthreads != -99) {
    bool newSetting = (nthreads != 1);
    if (!setting && newSetting) {
      setting = newSetting;
    }
  }
  return setting;
}

/// Small wrapper for tbb::task_arena.
/// Note that the tbbWrap::task_arena constructor is not thread-safe.
/// That should be fine because the task_arena is initialised before spawning
/// any threads.
class task_arena {
  std::optional<tbb::task_arena> tbb;

 public:
  task_arena(int nthreads = tbb::task_arena::automatic, unsigned res = 1) {
    if (enableTBB(nthreads)) {
      tbb.emplace(nthreads, res);
    }
  }

  template <typename F>
  void execute(const F& f) {
    if (tbb) {
      tbb->execute(f);
    } else {
      f();
    }
  }
};

/// Small wrapper for tbb::parallel_for.
class parallel_for {
 public:
  template <typename R, typename F>
  parallel_for(const R& r, const F& f) {
    if (enableTBB()) {
      tbb::parallel_for(r, f);
    } else {
      for (auto i = r.begin(); i != r.end(); ++i) {  // use default grainsize=1
        f(R(i, i + 1));
      }
    }
  }
};

/// Small wrapper for tbb::queuing_mutex and tbb::queuing_mutex::scoped_lock.
class queuing_mutex {
  std::optional<tbb::queuing_mutex> tbb;

 public:
  queuing_mutex() {
    if (enableTBB()) {
      tbb.emplace();
    }
  }

  class scoped_lock {
    std::optional<tbb::queuing_mutex::scoped_lock> tbb;

   public:
    scoped_lock() {
      if (enableTBB()) {
        tbb.emplace();
      }
    }

    scoped_lock(queuing_mutex& m) {
      if (enableTBB()) {
        tbb.emplace(*m.tbb);
      }
    }
  };
};

}  // namespace tbbWrap
