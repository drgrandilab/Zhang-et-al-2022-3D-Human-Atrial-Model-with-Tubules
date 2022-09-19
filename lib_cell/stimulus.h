#ifndef STIMULUS_H
#define STIMULUS_H

#include <math.h>
#include <stdio.h>
// #include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>

double S1(double stim_start, double stim, double BCL, double current_time, double stim_duration);

class StimFromInputFile
{
public:
    StimFromInputFile(const char *, bool report=false);
    ~StimFromInputFile() {};
    double ApplyStim(double stim_current, double stim_duration, double time_start, double current_time) ;
    std::vector<double> CL;
    double BCL;
    int StimNum, TotalNum;
    double StartTime, TimeSinceStim;
    bool stim_ON, stim_ON_prev;
    bool _inprogram_report;
};

#endif