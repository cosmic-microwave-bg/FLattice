#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <memory>
#include <vector>
#include <random>
#include "parameter.hpp"
#include "common.hpp"
#include "create_evolution_scheme.hpp"
#include "create_physical_quantity.hpp"



class Simulator
{
    private:
        Field _field;
        std::string _model_name;
        std::ofstream _status;
        std::unique_ptr<EvolutionScheme> _scheme;
        std::unordered_map<std::string,std::unique_ptr<PhysicalQuantity>> _quantity;
        std::vector<std::string> _quantity_name;

    public:
        Simulator( int num_fields, int N, int dimension )
        : _field(num_fields, N, dimension), _status(statusname.c_str(), std::ios::trunc), _scheme(nullptr)
        {
            if( !_status ){ std::cerr << "Failed to open 'status.txt'" << std::endl; exit(1); }
            else{
                _status << "    t";
                if( expansion != Expansion::no_expansion ) _status << "    a" << "    1/H";
                for( int n = 0; n < num_fields; ++n ) _status << "  field_ave["  << n << "]";
                for( int n = 0; n < num_fields; ++n ) _status << "  field_var["  << n << "]";
            }
        }

        void initializeField( int rnd ) {
            std::mt19937 mt(rnd);
            std::uniform_real_distribution<> rand(-1.e-5, 1.e-5);
    
            for( int n = 0; n < num_fields; ++n ){
                int i = 0, j = 0, k = 0;
                // Don't parallerize to realize the same initial condition
                for( i = 0; i < N; ++i ){
                    int idx = i;
                    #if DIMENSION >= 2
                        for( j = 0; j < N; ++j ){
                            idx = i*N+j;
                    #endif
                    #if DIMENSION == 3
                            for( k = 0; k < N; ++k ){
                                idx = (i*N+j)*N+k;
                    #endif
                                double f_fluct = rand(mt);
                                double v_fluct = rand(mt);
                                _field.f[n][idx]   = ini_amp[n]*(1 + f_fluct);
                                _field.df[n][idx]  = ini_vel[n];
                    #if DIMENSION == 3
                            }
                    #endif
                    #if DIMENSION >= 2
                        }
                    #endif
                }
            }
        }

        bool setModel (std::string name) { _model_name = name; }

        bool setEvolutionScheme (std::string name, int precision, Expansion expansion)
        {
            _scheme = createEvolutionScheme(name, precision, expansion, _field.f, _field.df);
            if( !_scheme ){ std::cerr << "invalid evolution scheme" << std::endl; exit(1); }
        }
        
        bool addPhysicalQuantites (std::string name)
        {
            _quantity[name] = createPhysicalQuantity(name, _field.f, _field.df);
            if( _quantity[name] == nullptr ){ std::cerr << "invalid physical quantity " << std::endl; exit(1); }
            _quantity_name.push_back(name);
            for( int n = 0; n < num_fields; ++n ) _status << "  "+name+"_ave[" << n << "] ";
            _status << "  total_"+name+"_ave";
        }
        
        void run( int output_step ) { _scheme->evolution(_field.f, _field.df, output_step); }

        void calculatePhysicalQuantities()
        { for( auto itr = _quantity.begin(); itr != _quantity.end(); ++itr ) itr->second->calculate(_field.f, _field.df); }

        void writeFields( int loop=0 )
        { writeVTI<double>( _field.f, num_fields, _model_name, loop ); }

        void writePhysicalQuantities( int loop=0, std::string name="none", bool write_each=true, bool write_total=true )
        {
            for( auto itr = _quantity.begin(); itr != _quantity.end(); ++itr ){
                std::string qname = itr->first;
                if( name == "none" or name == qname ){
                    if( write_each ) writeVTI<double>( itr->second->data(), num_fields, qname, loop );
                    if( write_total ) writeVTI<double>( itr->second->data_tot(), qname+"_total", loop );
                }
            }
        }

        void writeStatus()
        {
            double **f  = _field.f;
            double **df = _field.df;

            _status << std::endl;

            std::vector<double> field_ave(num_fields);
            std::vector<double> field_var(num_fields);
            for( int n = 0; n < num_fields; ++n ){
                field_ave[n] = calculateAverage(f[n], N, DIMENSION);
                field_var[n] = calculateVariance(f[n], field_ave[n], N, DIMENSION);
            }
            
            _status << std::noshowpos << std::fixed <<std::setprecision(2);
            _status << t;
            
            if( expansion != Expansion::no_expansion ){
                _status << "  " << a;
                if(expansion == Expansion::rad_domination) _status << "  " << ((A==0)? 2*t: c*pow(t, 2));
                if(expansion == Expansion::mat_domination) _status << "  " << ((A==0)? 3./2*t: c*pow(t/2, 3));
                if(expansion == Expansion::self_consist) _status << "  " << 1/(sqrt(calculateAverage(_quantity["energy"]->data_tot(), N, DIMENSION)*2/(D*(D-1))) * R);
            }

            _status << std::showpos << std::scientific << std::setprecision(3);
            for( int n = 0; n < num_fields; ++n ) _status << "  " << field_ave[n];
            for( int n = 0; n < num_fields; ++n ) _status << "  " << field_var[n];

            for( auto name: _quantity_name ){
                for( int n = 0; n < num_fields; ++n )
                    _status << "  " << calculateAverage(_quantity[name]->data(n), N, DIMENSION);
                _status << "  " << calculateAverage(_quantity[name]->data_tot(), N, DIMENSION);
            }
        }

        template <typename T=double>
        void searchObject( double threshold, int loop )
        {
            std::stringstream filename1;
            filename1 << "../data/QEsize_th=";
            filename1 << std::scientific << std::setprecision(0) << threshold << "_" << std::setw(4) << std::setfill('0') << loop << ".txt";

            auto QEs = serachObject<T>( _quantity["charge"]->data_tot(), _quantity["energy"]->data_tot(), threshold );
            std::ofstream ofs( filename1.str().c_str(), std::ios::app );
            for( auto p: QEs ){
                ofs << rnd << " ";
                for( auto v: p ) ofs << v << " ";
                ofs << std::endl;
            }
            ofs << std::endl;
        }

};



#endif