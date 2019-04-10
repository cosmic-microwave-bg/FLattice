#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include <chrono>


class Stopwatch
{
    private:
    	std::chrono::high_resolution_clock::time_point start, rec;
    
    public:
        Stopwatch():
        start(std::chrono::high_resolution_clock::now()), rec(std::chrono::high_resolution_clock::now()) {}

        void   set() { start = rec = std::chrono::high_resolution_clock::now(); }
        double lap()
        {
            if( rec == start ){
                rec = std::chrono::high_resolution_clock::now();
                return std::chrono::duration_cast<std::chrono::milliseconds>(rec - start).count()*1.e-3;
            }else{

                std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
                double tmp = std::chrono::duration_cast<std::chrono::milliseconds>(now - rec).count()*1.e-3;
                rec = now;
                return tmp;
            }
        }
        double end()
        {
            rec = std::chrono::high_resolution_clock::now();
            return std::chrono::duration_cast<std::chrono::milliseconds>(rec - start).count()*1.e-3;
        }
};



#endif