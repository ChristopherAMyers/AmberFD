#include "timer.h"

Timer::Timer()
{}
Timer::~Timer()
{}

void Timer::add_timer(TimerNames timer_name)
{
    timers.insert(std::pair<TimerNames, TimerInfo>(timer_name, TimerInfo()));
}

void Timer::start(TimerNames timer_name)
{
    if(timers.find(timer_name) != timers.end())
        add_timer(timer_name);
    timers[timer_name].start = omp_get_wtime();

}
void Timer::stop(TimerNames timer_name)
{
    timers[timer_name].total += omp_get_wtime() - timers[timer_name].start;
    timers[timer_name].n_calls += 1;
}

void Timer::print_timings()
{
    printf("\n\n TIMINGS INFO\n");
    // std::map<TimerNames, double>::iterator it;
    // for(it = timers.begin(); it != timers.end(); it++)
    for (auto const& it: timers)
    {
        TimerNames timer_name = it.first;
        double total = it.second.total;
        //printf("     %20s: %15.10f \n", timer_name.c_str(), total);
    }
}