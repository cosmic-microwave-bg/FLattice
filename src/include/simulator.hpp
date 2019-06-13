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
#include "create_model.hpp"



class Simulator
{
    private:
        Field _field;
        std::ofstream _status;
        std::shared_ptr<Model> _model;
        std::unique_ptr<EvolutionScheme> _scheme;
        std::unordered_map<std::string,std::unique_ptr<PhysicalQuantity>> _quantity;
        std::vector<std::string> _quantity_name;

    public:
        Simulator( int num_fields, int N, int dimension )
        : _field(num_fields, N, dimension), _status("../status.txt", std::ios::trunc), _model(nullptr), _scheme(nullptr)
        {
            if( !_status ){ std::cerr << "Failed to open 'status.txt'" << std::endl; exit(1); }
            else{
                _status << "    t";
                if( expansion != Expansion::no_expansion ) _status << "    a";
                for( int n = 0; n < num_fields; ++n ) _status << "  field_ave["  << n << "]";
                for( int n = 0; n < num_fields; ++n ) _status << "  field_var["  << n << "]";
            }
        }

        void initializeField( int rnd ) {
            double **f  = _field.f;
            double **df = _field.df;

            std::mt19937 mt( rnd );
            std::uniform_real_distribution<> rand(-1.e-5, 1.e-5);
    
            for( int n = 0; n < num_fields; ++n ){
                #pragma omp parallel for schedule( static ) num_threads ( num_threads )
                for( int i = 0; i < N; ++i ){	
                    #pragma omp simd
                    for( int j = 0; j < N; ++j ){
                        int idx = i*N+j;
                        double f_fluct = rand(mt);
                        double v_fluct = rand(mt);
                        f[n][idx]  = 10*(1 + f_fluct);
                        if( n == 0 ) df[n][idx] = f[n][idx];
                        if( n == 1 ) df[n][idx] = f[n][idx] + 1*(1 + v_fluct);
                    }
                }
            }
        }

        bool setModel (std::string name)
        {
            _model = createModel(name);
            if( !_model ){
                std::cerr << "invalid model name: " << name << std::endl;
                exit(1);
            }
        }

        bool setEvolutionScheme (std::string name, int precision, Expansion expansion)
        {
            _scheme = createEvolutionScheme(name, precision, _model, expansion, _field.f, _field.df);
            if( !_scheme ){ std::cerr << "invalid evolution scheme" << std::endl; exit(1); }
        }
        
        bool addPhysicalQuantites (std::string name)
        {
            _quantity[name] = createPhysicalQuantity(name, _model, _field.f, _field.df);
            if( _quantity[name] == nullptr ){ std::cerr << "invalid physical quantity " << std::endl; exit(1); }
            _quantity_name.push_back(name);
            for( int n = 0; n < num_fields; ++n ) _status << "  "+name+"_ave[" << n << "] ";
            _status << "  total_"+name+"_ave";
        }
        
        void run( int output_step ) { _scheme->evolution(_field.f, _field.df, output_step); }

        void calculatePhysicalQuantities()
        { for( auto itr = _quantity.begin(); itr != _quantity.end(); ++itr ) itr->second->calculate(_field.f, _field.df); }

        void writeFields( int loop=0 )
        { writeVTI<double>( _field.f, num_fields, _model->name(), loop ); }

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
                for( int n = 0; n < num_fields; ++n ) field_ave[n] /= a;
                for( int n = 0; n < num_fields; ++n ) field_var[n] /= a*a;
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
};



#endif