//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//


#ifndef THREAD_UTILS
#define THREAD_UTILS

#include <thread>
#include <queue>
#include <future>
#include <functional>
#include <atomic>

#define USE_TBB

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

namespace threadutils {
  inline unsigned get_num_threads()
  {
#ifdef NDEBUG
    return std::thread::hardware_concurrency();
#else
    return 1U;
#endif
  }
  
  template<typename T>
  class thread_safe_queue
  {
  private:
    mutable std::mutex mut;
    std::queue<std::shared_ptr<T> > data_queue;
    std::condition_variable data_cond;
  public:
    thread_safe_queue()
    {}
    
    void wait_and_pop(T& value)
    {
      std::unique_lock<std::mutex> lk(mut);
      data_cond.wait(lk,[this]{return !data_queue.empty();});
      value=std::move(*data_queue.front());  // 1
      data_queue.pop();
    }
    
    bool try_pop(T& value)
    {
      std::lock_guard<std::mutex> lk(mut);
      if(data_queue.empty())
        return false;
      value=std::move(*data_queue.front());  // 2
      data_queue.pop();
      return true;
    }
    
    std::shared_ptr<T> wait_and_pop()
    {
      std::unique_lock<std::mutex> lk(mut);
      data_cond.wait(lk,[this]{return !data_queue.empty();});
      std::shared_ptr<T> res=data_queue.front();  // 3
      data_queue.pop();
      return res;
    }
    
    std::shared_ptr<T> try_pop()
    {
      std::lock_guard<std::mutex> lk(mut);
      if(data_queue.empty())
        return std::shared_ptr<T>();
      std::shared_ptr<T> res=data_queue.front();  // 4
      data_queue.pop();
      return res;
    }
    
    void push(T new_value)
    {
      std::shared_ptr<T> data(
                              std::make_shared<T>(std::move(new_value)));  // 5
      std::lock_guard<std::mutex> lk(mut);
      data_queue.push(data);
      data_cond.notify_one();
    }
    
    bool empty() const
    {
      std::lock_guard<std::mutex> lk(mut);
      return data_queue.empty();
    }
  };
  
  class join_threads
  {
    std::vector<std::thread>& threads;
  public:
    explicit join_threads(std::vector<std::thread>& threads_):
    threads(threads_)
    {}
    ~join_threads()
    {
      for(unsigned long i=0;i<threads.size();++i)
      {
        if(threads[i].joinable())
          threads[i].join();
      }
    }
  };
  
  class function_wrapper
  {
    struct impl_base {
      virtual void call()=0;
      virtual ~impl_base() {}
    };
    
    std::unique_ptr<impl_base> impl;
    template<typename F>
    struct impl_type: impl_base
    {
      F f;
      impl_type(F&& f_): f(std::move(f_)) {}
      void call() { f(); }
    };
  public:
    template<typename F>
    function_wrapper(F&& f):
    impl(new impl_type<F>(std::move(f)))
    {}
    
    void operator()() { impl->call(); }
    
    function_wrapper() = default;
    
    function_wrapper(function_wrapper&& other):
    impl(std::move(other.impl))
    {}
    
    function_wrapper& operator=(function_wrapper&& other)
    {
      impl=std::move(other.impl);
      return *this;
    }
    
    function_wrapper(const function_wrapper&)=delete;
    function_wrapper(function_wrapper&)=delete;
    function_wrapper& operator=(const function_wrapper&)=delete;
  };
  
  class thread_pool
  {
    std::atomic_bool done;
    thread_safe_queue<function_wrapper> work_queue;  // 使用function_wrapper，而非使用std::function
    std::vector<std::thread> threads;  // 2
    join_threads joiner;  // 3
    
    void worker_thread()
    {
      while(!done)
      {
        function_wrapper task;
        if(work_queue.try_pop(task))
        {
          task();
        }
        else
        {
          std::this_thread::yield();
        }
      }
    }
  public:
    const std::vector<std::thread>& get_threads() const
    {
      return threads;
    }
    
    thread_pool():
    done(false),joiner(threads)
    {
      unsigned const thread_count = get_num_threads();
      try
      {
        for(unsigned i=0;i<thread_count;++i)
        {
          threads.push_back(std::thread(&thread_pool::worker_thread,this));  // 9
        }
      }
      catch(...)
      {
        done=true;  // 10
        throw;
      }
    }
    
    ~thread_pool()
    {
      done=true;  // 11
    }
    
    template<typename FunctionType>
    std::future<typename std::result_of<FunctionType()>::type>  // 1
    submit(FunctionType f)
    {
      typedef typename std::result_of<FunctionType()>::type
      result_type;  // 2
      
      std::packaged_task<result_type()> task(std::move(f));  // 3
      std::future<result_type> res(task.get_future());  // 4
      work_queue.push(std::move(task));  // 5
      return res;  // 6
    }
    
    template<typename Index, typename Callable>
    static void ParallelFor(Index start, Index end, Callable func) {
#ifdef NDEBUG
#ifdef USE_TBB
      tbb::parallel_for(start, end, 1, func);
#else
      // Estimate number of threads in the pool
      const static unsigned nb_threads_hint = std::thread::hardware_concurrency();
      const static unsigned nb_threads = (nb_threads_hint == 0u ? 8u : nb_threads_hint);
      
      // Size of a slice for the range functions
      Index n = end - start + 1;
      Index slice = (Index) std::round(n / static_cast<double> (nb_threads));
      slice = std::max(slice, Index(1));
      
      // [Helper] Inner loop
      auto launchRange = [&func] (int k1, int k2) {
        for (Index k = k1; k < k2; k++) {
          func(k);
        }
      };
      
      // Create pool and launch jobs
      std::vector<std::thread> pool;
      pool.reserve(nb_threads);
      Index i1 = start;
      Index i2 = std::min(start + slice, end);
      for (unsigned i = 0; i + 1 < nb_threads && i1 < end; ++i) {
        pool.emplace_back(launchRange, i1, i2);
        i1 = i2;
        i2 = std::min(i2 + slice, end);
      }
      if (i1 < end) {
        pool.emplace_back(launchRange, i1, end);
      }
      
      // Wait for jobs to finish
      for (std::thread &t : pool) {
        if (t.joinable()) {
          t.join();
        }
      }
#endif
#else
      for(Index i = start; i < end; ++i) {
        func(i);
      }
#endif
    }
  };
};
#endif /* ThreadUtils_hpp */
