#ifndef _SEARCH_H_
#define _SEARCH_H_

#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <memory>
#include <vector>
#include <queue>
#include "parameter.hpp"
#include "common.hpp"


// chargeに対しthreshold以上の領域のchageとenergyを計算し、そのvectorを返す
// Only 2D for now
template <typename T=double>
std::vector<std::vector<T>> serachObject( double* charge, double* energy, double threshold )
{
    std::vector<bool> is_visited(pow(N,DIMENSION), false);

    std::vector<std::vector<T>> qball;
    for( int i = 0; i < N; ++i ){
        for( int j = 0; j < N; ++j ){
            int idx = i*N+j;
            if( is_visited[idx] or charge[idx] < threshold ) continue;

            double Q = 0, E = 0;
            int count = 0;
            std::queue<std::pair<int,int>> q;
            q.push( std::make_pair(i, j) );
            while( !q.empty() )
            {
                auto tmp = q.front(); q.pop();
                int x = tmp.first;
                int y = tmp.second;
                int now = x*N+y;

                if( !is_visited[now] and charge[now] >= threshold ){
                    is_visited[now] = true;
                    Q += charge[now];
                    E += energy[now];
                    ++count;

                    int xp1 = (x == N-1)?     0: x+1;
                    int xm1 = (x ==   0)?   N-1: x-1;
                    int yp1 = (y == N-1)?     0: y+1;
                    int ym1 = (y ==   0)?   N-1: y-1;

                    int nx[8] = {  x, xp1, xp1, xp1,   x, xm1, xm1, xm1};
                    int ny[8] = {yp1, yp1,   y, ym1, ym1, ym1,   y, yp1};
                    for( int l = 0; l < 8; ++l ){
                        int next = nx[l]*N+ny[l];
                        if( !is_visited[next] and charge[next] >= threshold )
                            q.push( std::make_pair(nx[l], ny[l]) );
                    }
                }
            }
            qball.push_back( {Q*pow(a*dx,D), E*pow(a*dx,D), (double)count} );  // [charge, energy, size], for physical size: sqrt(count/(4*M_PI))*dx
        }
    }
    return qball;
}



#endif