#include "threadpool.hpp"

#include <stdexcept>
#include <utils/logger.hpp>

namespace utils {

ThreadPool::ThreadPool(size_t threads_num, const std::string& name)
    : threads_num_(threads_num) {
  Run();
}

ThreadPool::~ThreadPool() { Stop(); }

std::future<void> ThreadPool::Execute(ThreadPool::Task task) {
  auto wrapper =
      std::make_shared<std::packaged_task<decltype(task())()>>(std::move(task));
  {
    std::unique_lock<std::mutex> lock(eventMutex_);
    if (stopping_) {
      throw std::runtime_error("ThreadPool has been stopped");
    }
    tasks_.emplace([=] { (*wrapper)(); });
  }

  event_.notify_one();
  return wrapper->get_future();
}

size_t ThreadPool::GetQueueSize() const {
  std::unique_lock<std::mutex> lock(eventMutex_);
  return tasks_.size();
}

size_t ThreadPool::GetThreadCount() const { return threads_.size(); }

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
    threads_.emplace_back([this] {
      while (true) {
        Task task;
        {
          std::unique_lock<std::mutex> lock{eventMutex_};
          event_.wait(lock, [this] { return stopping_ || !tasks_.empty(); });

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
