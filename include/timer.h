#include <stdio.h>
#include <map>
#include <string>
#include <omp.h>

#ifndef TIMER
#define TIMER
class TimerInfo
{
    public:
        TimerInfo(): total(0.0), n_calls(0)
        {}
        double start;
        double total;
        unsigned int n_calls;
};

enum TimerNames{
    threaded_energy,
    energy_pairs,
    minimization,
    dac_fragments,
    dac_total
};

class Timer
{
    public:
        Timer();
        ~Timer();
        void add_timer(TimerNames timer_name);
        void start(TimerNames timer_name);
        void stop(TimerNames timer_name);
        void print_timings();
    private:
        std::map<TimerNames, TimerInfo> timers;
};
#endif // TIMER