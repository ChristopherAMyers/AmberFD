#include <stdio.h>
#include <map>
#include <string>
#include <omp.h>

#ifndef TIMER_INFO
#define TIMER_INFO
class TimerInfo
{
    public:
        TimerInfo(): total(0.0), n_calls(0)
        {}
        double start;
        double total;
        unsigned int n_calls;
};
#endif // TIMER_INFO

#ifndef TIMER
#define TIMER
class Timer
{
    public:
        Timer();
        ~Timer();
        void add_timer(std::string timer_name);
        void start(std::string timer_name);
        void stop(std::string timer_name);
        void print_timings();
    private:
        std::map<std::string, TimerInfo> timers;
};
#endif // TIMER