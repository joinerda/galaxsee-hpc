#ifndef TIMER_H
#define TIMER_H

#include <unistd.h>
#include <time.h>

typedef struct {
    struct timespec begin;
    struct timespec end;
    double elapsed_time;
    int counter;
} Timer;

void TimerReset(Timer * timer);
void TimerStart(Timer * timer);
void TimerStop(Timer * timer);
double TimerRead(Timer * timer);
double TimerReadPer(Timer * timer);
int TimerCount(Timer * timer);

#endif
