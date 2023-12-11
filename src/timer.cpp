#include "timer.h"
#include <chrono>
#include <vector>
using namespace std; 

size_t timer::n_now = 0;
std::map<std::string,std::map<std::string,timer::MyTimer>> timer::timer_pool;
timer::timer(){} 
timer::~timer(){}

void timer::finish(std::ofstream &ofs,const bool print_flag)
{
        timer::tick("","total");
        if(print_flag)
                timer::print_all(ofs);
}

void timer::start()
{
        timer::tick("","total");
        return ;
}

double timer::cpu_time(void)
{
	//t2-t1
        static auto t1 = std::chrono::system_clock::now();
        const auto t2 = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        return double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
}

void timer::tick(const std::string& class_name,const std::string& function_name)
{
        	MyTimer &mytimer = timer_pool[class_name][function_name];
		if(mytimer.start_flag){
			mytimer.cpu_start = cpu_time();
			++mytimer.calls;
			mytimer.start_flag = false;
		}	
		else{
			mytimer.cpu_second += (cpu_time() - mytimer.cpu_start);
			mytimer.start_flag = true;
		}
}

void timer::print_all(std::ofstream &ofs)
{      
        std::vector<std::pair<std::pair<std::string,std::string>,MyTimer>> timer_pool_order;
        for(auto &timer_pool_A : timer_pool)
        {
                const std::string class_name = timer_pool_A.first;
                for(auto &timer_pool_B : timer_pool_A.second)
                {
                        const std::string name = timer_pool_B.first;
                        const MyTimer mytimer = timer_pool_B.second;
                        if(timer_pool_order.size() < mytimer.order+1)
                                timer_pool_order.resize(mytimer.order+1);
                        timer_pool_order[mytimer.order] = std::pair<std::pair<std::string,std::string>, MyTimer> {std::pair<std::string,std::string >{class_name,name}, mytimer};
                }
        }

        std::cout << std::setprecision(2);
        ofs << std::setprecision(3);
        std::cout<<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG--------|PER%-------" << std::endl;
        ofs <<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG--------|PER%-------" << std::endl;
        for(auto &timer_pool_order_A : timer_pool_order)
        {
                const std::string &class_name = timer_pool_order_A.first.first;
                const std::string &function_name = timer_pool_order_A.first.second;
                const MyTimer &mytimer = timer_pool_order_A.second;
                
                ofs << std::resetiosflags(std::ios::scientific);
                ofs  << " "
                         << std::setw(2)  << " "
                         << std::setw(20) << class_name
                         << std::setw(20) << function_name
                         << std::setw(15) << std::setprecision(10) << mytimer.cpu_second
                         << std::setw(10) << mytimer.calls
                         << std::setw(10) << std::setprecision(10) << mytimer.cpu_second/mytimer.calls
                         << std::setw(20) << mytimer.cpu_second / timer_pool_order[0].second.cpu_second * 100 << "%" << std::endl;
                
                std::cout << std::resetiosflags(std::ios::scientific);

                std::cout << " "
                         << std::setw(2)  << " "
                         << std::setw(20) << class_name
                         << std::setw(20) << function_name
                         << std::setw(15) << std::setprecision(6) << mytimer.cpu_second
                         << std::setw(10) << mytimer.calls
                         << std::setw(10) << std::setprecision(6) << mytimer.cpu_second/mytimer.calls
                         << std::setw(15) << mytimer.cpu_second / timer_pool_order[0].second.cpu_second * 100 << "%" << std::endl;
        }
        std::cout<<" ----------------------------------------------------------------------------------------"<<std::endl;
        ofs <<" ----------------------------------------------------------------------------------------"<<std::endl;
}
 






