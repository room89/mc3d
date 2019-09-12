#include "threadpool.hpp"

#include <utils/logger.hpp>

namespace utils {

ThreadPool::ThreadPool(size_t threads_num, const std::string& name)
    : threads_num_(threads_num) {
  Run();
}

ThreadPool::~ThreadPool() {}

std::future<void> ThreadPool::Execute(ThreadPool::Task task) {
  auto wrapper =
      std::make_shared<std::packaged_task<decltype(task())()>>(std::move(task));
  {
    std::unique_lock<std::mutex> lock(eventMutex_);
    tasks_.emplace([=] { (*wrapper)(); });
  }

  event_.notify_one();
  return wrapper->get_future();
}

size_t ThreadPool::GetQueueSize() const { return tasks_.size(); }

void ThreadPool::Run() {
  threads_.reserve(threads_num_);
  for (size_t i = 0; i < threads_num_; ++i) {
    threads_.emplace_back([this] {
      while (true) {
        Task task;
        {
          std::unique_lock<std::mutex> lock{eventMutex_};
          event_.wait(lock, [=] { return stopping_ || !tasks_.empty(); });

          if (stopping_) break;

          task = std::move(tasks_.front());
          tasks_.pop();
        }
        task();
      }
    });
  }
}

}  // namespace utils
