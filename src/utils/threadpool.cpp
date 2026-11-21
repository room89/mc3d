#include "threadpool.hpp"

#include <stdexcept>
#include <utils/logger.hpp>

namespace utils {

ThreadPool::ThreadPool(size_t threads_num, const std::string& name)
    : threads_num_(threads_num), active_threads_(threads_num) {
  Run();
}

ThreadPool::~ThreadPool() { Stop(); }

std::future<void> ThreadPool::Execute(ThreadPool::Task task) {
  auto wrapper =
      std::make_shared<std::packaged_task<decltype(task())()>>(std::move(task));
  bool limited;
  {
    std::unique_lock<std::mutex> lock(eventMutex_);
    if (stopping_) {
      throw std::runtime_error("ThreadPool has been stopped");
    }
    tasks_.emplace([=] { (*wrapper)(); });
    limited = active_threads_ != threads_num_;
  }

  if (limited) event_.notify_all();
  else event_.notify_one();
  return wrapper->get_future();
}

size_t ThreadPool::GetQueueSize() const {
  std::unique_lock<std::mutex> lock(eventMutex_);
  return tasks_.size();
}

size_t ThreadPool::GetThreadCount() const {
  std::lock_guard<std::mutex> lock(eventMutex_);
  return active_threads_;
}

void ThreadPool::SetActiveThreadCount(size_t count) {
  {
    std::lock_guard<std::mutex> lock(eventMutex_);
    if (count == 0 || count > threads_num_)
      throw std::invalid_argument("Active thread count outside pool capacity");
    if (active_threads_ == count) return;
    active_threads_ = count;
  }
  event_.notify_all();
}

void ThreadPool::Stop() {
  {
    std::unique_lock<std::mutex> lock(eventMutex_);
    if (stopping_) {
      return;
    }
    stopping_ = true;
  }
  event_.notify_all();
  for (auto& thread : threads_) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void ThreadPool::Run() {
  threads_.reserve(threads_num_);
  for (size_t i = 0; i < threads_num_; ++i) {
    threads_.emplace_back([this, i] {
      while (true) {
        Task task;
        {
          std::unique_lock<std::mutex> lock{eventMutex_};
          event_.wait(lock, [this, i] {
            return stopping_ || (i < active_threads_ && !tasks_.empty());
          });

          if (stopping_ && tasks_.empty()) break;

          task = std::move(tasks_.front());
          tasks_.pop();
        }
        task();
      }
    });
  }
}

}  // namespace utils
