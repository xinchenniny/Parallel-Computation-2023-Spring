#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

using namespace std;

class  timer{
public:
	timer();
	~timer();
    struct MyTimer  
    {
        double cpu_start;
        double cpu_second = 0.0;
        size_t calls = 0;
        size_t order = n_now++;
        bool start_flag = true;
    };
	static std::map<std::string, std::map<std::string, MyTimer>> timer_pool;

	static void tick (const std::string& class_name,const std::string& function_name);
	
	static void start(void);

	static void finish(std::ofstream &ofs, const bool print_flag);
	// to control whether to use tick 
	static void on(void)
    {
	  flag_for_tick = false;
    }
    
    static void off(void)
    {
        flag_for_tick= true;
    }
    
    static void print_all(std::ofstream &ofs);
    
    static long double print_until_now(void);

private:
    static bool flag_for_tick;
    static size_t n_now;
    static double cpu_time(void);

};
#endif
