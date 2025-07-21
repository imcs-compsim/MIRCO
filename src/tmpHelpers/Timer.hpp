#pragma once

#include <chrono>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#define TIMERSON true

class Timer;  // fwd decl

// GLOBAL DEFINITIONS {

extern int numNonlinearIters;
extern int numLinearCallsTotal;
extern int linearSystemSizeSummed;
extern double linearSystemSizeAvg;
extern int evaluateIterations;

//}

#if (TIMERSON)

class TimerRegistry
{
 public:
  /// TimerRegistry() {};
  TimerRegistry(const std::string& registryName) : registryName_(registryName) {};

 public:
  static TimerRegistry& globalInstance()
  {
    static TimerRegistry registry = TimerRegistry("Global Instance");
    return registry;
  }

  // for global tracking
  void start() { totalStart_ = std::chrono::high_resolution_clock::now(); }

  // Usually called by the Timers themselves (automatically for ScopedTimer)
  void addTimer(const Timer& timer);

  std::string timingReportStr(bool sortedByFinishElseStart = false);

 private:
  std::string registryName_;
  mutable std::mutex mutex_;
  // Array<Timer> timers_;
  std::chrono::high_resolution_clock::time_point totalStart_;
  std::unordered_map<std::string, Timer> timers_;
  std::vector<std::string> timersOrdered_;
};



class Timer
{  // TODO:m probably should all virtual destructor for memory stuff? so the destruction gets
   // upcast? idk how it works honestly
 public:
  Timer() {}
  Timer(const std::string& name) : name_(name) {}
  double getTimerVal_secDouble() const { return timerVal_secDouble_; }

 public:
  std::string name_;
  std::chrono::high_resolution_clock::time_point start_;
  double timerVal_secDouble_ = 0.0;

 protected:
  bool active_ = false;
};

class ScopedTimer : public Timer
{
 public:
  ScopedTimer(const std::string& name) : Timer(name)
  {
    active_ = true;
    start_ = std::chrono::high_resolution_clock::now();
  }

  ~ScopedTimer()
  {
    active_ = false;
    timerVal_secDouble_ =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_).count();
    TimerRegistry::globalInstance().addTimer(*this);
  }
};

class StandardTimer : public Timer
{
 public:
  StandardTimer(const std::string& name) : Timer(name) { active_ = false; }

  void start()
  {
    active_ = true;
    start_ = std::chrono::high_resolution_clock::now();
    lastStart_ = start_;
  }
  void stop()
  {
    active_ = false;
    timerVal_secDouble_ +=
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - lastStart_)
            .count();
    TimerRegistry::globalInstance().addTimer(
        *this);  // TODO: this can be organized and done better with hybrid pause maybe for static
  }
  void pause()
  {
    timerVal_secDouble_ +=
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - lastStart_)
            .count();
  }
  void unpause() { lastStart_ = std::chrono::high_resolution_clock::now(); }

 private:
  std::chrono::high_resolution_clock::time_point lastStart_;
};

#else

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wunused-parameter"

class TimerRegistry
{
 public:
  TimerRegistry(const std::string& registryName) noexcept : registryName_(registryName) {};

 public:
  static TimerRegistry& globalInstance()
  {
    static TimerRegistry registry = TimerRegistry("Global Instance");
    return registry;
  }

  // for global tracking
  void start() { totalStart_ = std::chrono::high_resolution_clock::now(); }

  void addTimer(const Timer&) noexcept {}

  std::string timingReportStr() const noexcept;

 private:
  std::string registryName_;
  std::chrono::high_resolution_clock::time_point totalStart_;
};

class Timer
{
 public:
  Timer() noexcept {}
  Timer(const std::string&) noexcept {}
  ~Timer() noexcept = default;
};

class ScopedTimer
{
 public:
  ScopedTimer(const std::string&) noexcept {}
  ~ScopedTimer() noexcept = default;
};


class StandardTimer
{
 public:
  StandardTimer(const std::string&) noexcept {}
  ~StandardTimer() noexcept = default;
};

#endif
