#include <math.h>
#include <stdio.h>
#include "stimulus.h"

double S1(double stim_start_time, double stim_current, double BCL, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) >= 0 ? (current_time - stim_start_time) : -1.0;
    remain       =  fmod(time_elapsed, BCL);
    if (remain >= 0 && remain < stim_duration) {
        Istim   = stim_current;
    } else {
        Istim   = 0.0;
    }

    return Istim;
}