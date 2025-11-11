#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace utils {

class ThreadPool {
 public:
  ThreadPool(size_t threads_num, const std::string& name);
  ~ThreadPool();

  using Task = std::function<void()>;
  using Future = std::future<void>;

  void Stop();

  std::future<void> Execute(Task task);

  size_t GetQueueSize() const;

  size_t GetThreadCount() const;

 private:
  void Run();

  const size_t threads_num_;
  std::vector<std::thread> threads_;
  std::queue<Task> tasks_;

  std::condition_variable event_;
  mutable std::mutex eventMutex_;
  bool stopping_ = false;
};

}  // namespace utils
