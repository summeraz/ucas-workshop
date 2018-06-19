/*
 *  relax.cpp
 *
 *  Created by Andrew Z. Summers on 6/18/18.
 *  Modified from code created by Christopher Iacovella on 8/1/11.
 */

#include "main.h"
#include "energy.hpp"
#include "io.hpp"
#include "neighbor.hpp"


void run(double sigma, double epsilon, double cutoff, const char* coord_file,
         double temperature, double dx, int n_mc, int output_freq)
{
	
	// System parameters
	double beta = 1.0/(kb*temperature);	
	
	/*  
        Define arrays to hold position (x)
	    Note `coordlist_t` is just a convenient container class I defined
        that is essentially a double array x[n_particles][3]
    */
	coordlist_t x;

	// Define our box array
	double L[3];

	// Load particle positions and box information
	load_gro(coord_file, x, L);	
		
    int n_particles = x.size();
    double number_density = n_particles / (L[0] * L[1] * L[2]);
	
	double pe[n_particles];
	
	/*
        Initialize the neighborlist routine
        Note, for a small systems sizes, a simple n^2 neighborlist is most
        efficient however, for larger systems, a cell list would provide
        higher performance.
        Also note, this is a "full" neighborlist, as opposed to the "half"
        neighborlist typically used in MD. This is necessary as each MC
        displacement is essentially an independent event.
    */
	std::vector<neighbor> nbr;
	double skin = 0.5;
	nsq_neighbor_init(x, nbr, cutoff, skin, L);
	
	// Create output streams to write out trajectory data and thermodynamic data
	std::ofstream dataOut("run.xyz");	
	std::ofstream thermoOut("run.log");	
    print_xyz(dataOut, x, L);

	// Number of accepted moves
	int n_accept=0; 

	// System potential
	double potential =0;	

    /*
        ===============================================
        Main production loop
        ===============================================
    */
       
    for(int time = 0; time < n_mc; time++)
    {   
        int n_accept_local;
       
        n_accept_local = pair_LJ_neighborlist(x, nbr, sigma, epsilon, cutoff, L,
                                              dx, &potential, beta);
        n_accept += n_accept_local;

        // Check to see if we need to update our neighborlist
        nsq_neighbor_check(x, nbr, cutoff, skin, L);

        // Check to see if we should write data to file and to the screen
        if(time % output_freq==0)
        {
            double prob = (double)n_accept / (double)(output_freq * n_particles);

            n_accept = 0;
            potential = 0;

            // Calculate the total system potential
            calc_pe(x, nbr, sigma, epsilon, cutoff, L, &potential);

            std::cout << time << "\tPE: " << potential / (double)n_particles <<  "\tacceptance probability: " << prob << std::endl;

            thermoOut << time << "\t" << potential / (double)n_particles << std::endl;
            // Write trajectory data
            print_xyz(dataOut, x, L);
        }
    }
}


int main(int argc, char **argv)
{
    // LJ Parameters
    double sigma = atof(argv[1]);
    double epsilon = atof(argv[2]);
    double cutoff = atof(argv[3]);

    // MC Parameters
	double temperature = atof(argv[5]);
	double dx = atof(argv[6]);
    int n_mc = atoi(argv[7]);
    int output_freq = atoi(argv[8]);

	run(sigma, epsilon, cutoff, "system_relaxed.gro", temperature, dx, n_mc,
        output_freq);
	
	return 0;
}
