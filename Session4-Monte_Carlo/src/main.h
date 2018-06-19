/*
 *  main.h
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Copyright 2011. All rights reserved.
 *
 */


#include <stdlib.h>
#include <stdio.h>

#include <cmath>
#include <vector>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>

#define kb 1

#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))

#define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)

//define a custom vector type to make managing particle coordinates easier
typedef std::vector<double> coord_t;
typedef std::vector<coord_t> coordlist_t;


//define stream operators for the coord_t and coordlist class to make outputting data easier
std::ostream& operator << (std::ostream&, coordlist_t&);
std::ostream& operator << (std::ostream&, coord_t&);

std::istream& operator >> (std::istream&, coordlist_t&);
std::istream& operator >> (std::istream&, coord_t&);

void print(coordlist_t&, std::ostream&);
void print(coord_t&, std::ostream&);


//function I/O prototypes

void load_raw(const char* filename, coordlist_t& x);
void print_xyz(const char* filename, coordlist_t& x, double *L);
void print_xyz(std::ostream& dataOut, coordlist_t& x, double *L);
void init_system(coordlist_t& x, int N, double density, double *L);


//neighborlist class holds a particle id and the initial position of a particle
class neighbor{

public:
	std::vector <int> member;
	coordlist_t x_old;
	
};



void nsq_neighbor_init(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L); 
void nsq_neighbor_check(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L);
void nsq_neighbor_rebuild(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L);

void calc_pe(coordlist_t& x, std::vector<neighbor> &nbr, double sigma, double epsilon, double cutoff, double *L, double *potential_energy_total);
int pair_LJ_neighborlist(coordlist_t& x, std::vector<neighbor> &nbr, double sigma, double epsilon, double cutoff, double *L, double dx, double *potential_energy_total, double _beta);
