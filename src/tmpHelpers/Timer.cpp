#include "Timer.hpp"


int numNonlinearIters = 0;
int numLinearCallsTotal = 0;  // should == numNonlinearIters?
int linearSystemSizeSummed = 0;
double linearSystemSizeAvg = 0.0;  // linearSystemSizeSummed/numLinearCallsTotal
int evaluateIterations = 0;

#if (TIMERSON)

#include <algorithm>

void TimerRegistry::addTimer(const Timer& timer)
{
  std::lock_guard<std::mutex> lock(mutex_);
  auto it = timers_.find(timer.name_);
  if (it == timers_.end())
  {
    timersOrdered_.push_back(timer.name_);
    timers_[timer.name_] = timer;  // this does copy
  }
  else
  {
    timers_[timer.name_].timerVal_secDouble_ +=
        timer.timerVal_secDouble_;  // direct use for efficiency?
  }
}



#define FOR(i, n) for (int i = 0; i < (n); ++i)


template <typename T>
class Array
{
  using size_t = std::size_t;

 private:
  size_t size_;
  std::vector<T> data_;

 public:
  // Default ctor
  Array() : size_(0), data_(0) {}
  // Standard constructor
  Array(size_t size) : size_(size), data_(size) {}
  // This constructor basically makes a zero vector with the same dimensions as other
  Array(const Array& other) : size_(other.size()), data_(other.raw()) {}
  // For literal
  Array(const std::vector<T>& raw) : size_(raw.size()), data_(raw) {}

  std::vector<T>& raw() { return data_; }
  const std::vector<T>& raw() const { return data_; }

  // Access by (i, j)
  T& operator()(size_t i) { return data_[i]; }
  const T& operator()(size_t i) const { return data_[i]; }

  size_t size() const { return size_; }

  // Find ele with the value and return the positions of all occurances; if not found, return empty
  // arr
  Array<size_t> find(T val) const
  {
    Array<size_t> arr = Array<size_t>();
    for (int i = 0; i < size_; ++i)
    {
      if ((*this)(i) == val) arr.push_back(i);
    }
    return arr;
  }

  // Expose some raw std::vector functions
  void resize(size_t newSize)
  {
    size_ = newSize;
    data_.resize(newSize);
  }
  /// void resize(size_t newSize, double val) {data_.resize(newSize);} //# I think unecessary
  void push_back(T ele)
  {
    size_++;
    data_.push_back(ele);
  };

  void deleteIndices(const Array<size_t>& indicesToDelete)
  {
    Array<size_t> keepIndices;
    for (size_t i = 0; i < size_; ++i)
    {
      if (indicesToDelete.find((size_t)i).size() == 0) keepIndices.push_back(i);
    }

    std::vector<T> newData(keepIndices.size());
    for (size_t iNew = 0; iNew < keepIndices.size(); ++iNew)
    {
      size_t iOld = keepIndices(iNew);
      newData[iNew] = (*this)(iOld);
    }

    size_ = keepIndices.size();
    data_ = std::move(newData);
  }
  void deleteIndices(const Array<int>& indicesToDelete)
  {  // TODO: think about which we want. size_t or int or both possibly but also maybe not very good
     // to have both
    Array<size_t> keepIndices;
    for (size_t i = 0; i < size_; ++i)
    {
      if (indicesToDelete.find((size_t)i).size() == 0) keepIndices.push_back(i);
    }

    std::vector<T> newData(keepIndices.size());
    for (size_t iNew = 0; iNew < keepIndices.size(); ++iNew)
    {
      size_t iOld = keepIndices(iNew);
      newData[iNew] = (*this)(iOld);
    }

    size_ = keepIndices.size();
    data_ = std::move(newData);
  }

  // extends if newRowCount > size_; else return and throw no errors
  void extend(size_t newRowCount)
  {
    if (newRowCount <= size_) return;

    // inefficient
    for (int i = size_; i < newRowCount; ++i) (*this).push_back(0.0);
    size_ = newRowCount;
  }
};



// // STRING UTILITIES INTENDED FOR EVENTUAL OUTPUT (TODOi: put these in myUtils eventually, but
// probably need to put Array and such in myUtils as well, which I should do anyways)
inline std::string repeatStr(const std::string& s, int count)
{
  if (count <= 0) return "";
  std::string out = s;
  FOR(i, count - 1)
  out += s;
  return out;
}

inline std::string levelizeString(const std::string& s, int level)
{
  if (level == 0) return s;
  // else
  std::string lvlStr = repeatStr("\t", level);

  std::string tmpS = lvlStr;
  FOR(i, s.length())
  {
    tmpS += s[i];
    if (s[i] == '\n' && i != s.length() - 1)  // -1 because we don't want to turn the last "\n" into
                                              // "\n\t", or else the next thing is affected
      tmpS += lvlStr;
  }
  return tmpS;
}

inline Array<std::string> strToStrArray(const std::string& s)
{
  auto out = Array<std::string>();
  if (s == "") return out;
  out.push_back("");
  size_t sLength = s.length();
  FOR(i, sLength)
  {
    if (i != sLength - 1 && s[i] == '\n')
      out.push_back("");
    else if (s[i] != '\n')
      out(out.size() - 1) += s[i];
  }
  return out;
}
inline std::string strArrayToStr(const Array<std::string>& a)
{
  std::string out = "";
  FOR(i, a.size())
  out += a(i) + "\n";
  return out;
}

inline Array<int> checkForIn(const std::string& checkFor, const std::string& checkIn)
{  // TODOi: you have to move Array into myUtils so that you can also move this back into myUtils
  auto out = Array<int>();
  std::string checkInTemp = checkIn;
  if (checkFor.length() <= checkIn.length())
  {
    for (int i = 0; i <= checkIn.length() - checkFor.length(); ++i)
    {
      // checkInTemp = deleteInterval(checkInTemp, i, checkInTemp.Length);
      FOR(j, checkFor.length())
      {
        if (checkFor[j] != checkIn[i + j]) break;
        if (j == checkFor.length() - 1)
        {
          out.push_back(i);
        }
      }
    }
  }
  return out;
}

// aligns the string rows at certain keys, like e.g. "|"" will be vertically aligned according to
// count (TODOm: make an overload for Array<std::string> where each ele is a row)
inline std::string alignStringAt(const std::string& s, const std::string& alignerKey)
{
  auto a = strToStrArray(s);
  const int aSize = a.size();
  auto maxSpacing = Array<int>();
  Array<Array<int>> CFIs(aSize);  // store for efficiency
  FOR(i, aSize)
  {
    auto cfi = checkForIn(alignerKey, a(i));
    CFIs(i) = cfi;
    if (cfi.size() == 0) continue;
    maxSpacing.extend(cfi.size() - maxSpacing.size());
    if (cfi(0) > maxSpacing(0)) maxSpacing(0) = cfi(0);
    FOR(j, cfi.size() - 1)
    if (cfi(j + 1) - cfi(j) > maxSpacing(j + 1)) maxSpacing(j + 1) = cfi(j + 1) - cfi(j);
  }
  auto aOut = Array<std::string>(aSize);
  FOR(i, aSize)
  {
    auto& cfi = CFIs(i);
    int dist = 0;
    int k = 0;
    FOR(j, a(i).length())
    {
      if (k < cfi.size())
        if (j == cfi(k))
        {
          int diff = maxSpacing(k) - (cfi(k) - dist);
          aOut(i) += repeatStr(" ", diff);
          dist = cfi(k);
          ++k;
        }
      aOut(i) += a(i)[j];
    }
  }
  return strArrayToStr(aOut);
}



std::string TimerRegistry::timingReportStr(bool sortByStartElseFinish)
{
  std::string out = "==== Timing Report ====\n";
  out += " - Registry name: " + registryName_ + "\n";
  double totalTime =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart_)
          .count();
  out += " - Total registry time: " + std::to_string(totalTime) + "s\n";

  if (sortByStartElseFinish)
    std::sort(timersOrdered_.begin(), timersOrdered_.end(),
        [&](const std::string& a, const std::string& b)
        { return timers_.at(a).start_ < timers_.at(b).start_; });

  std::string out2 =
      "Timer name | Timer value [s] | Share of registry total [%]\n";  // | Start of timer [s]\n";
  FOR(i, timersOrdered_.size())
  {
    const std::string& name = timersOrdered_.at(i);
    const Timer& timer = timers_.at(name);
    double timerVal = timer.getTimerVal_secDouble();
    double shareOfTotal = timerVal / totalTime;
    out2 +=
        name + " | " + std::to_string(timerVal) + " | " + std::to_string(shareOfTotal * 100) +
        /*" | "+std::to_string(std::chrono::duration<double>(timer.start_-totalStart_).count())+*/
        "\n";
  }

  out2 = alignStringAt(out2, "|");

  return out + out2;
}

#else

std::string TimerRegistry::timingReportStr() const noexcept
{
  std::string out = "==== Timing Report : Timers are off (TIMERSON == false) ====\n";
  out += "Registry name: " + registryName_ + "\n";
  double totalTime =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - totalStart_)
          .count();
  out += "\tTotal registry time | " + std::to_string(totalTime) + "s\n";
  return out;
}

#endif
