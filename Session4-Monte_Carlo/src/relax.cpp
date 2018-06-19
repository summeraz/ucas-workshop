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


void relax(double sigma, double epsilon, double cutoff, const char* coord_file,
           double temperature, double dx, double target, int seed, int n_relax,
           int adjust_freq)
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
	
	// Print out system parameters
	std::cout << "=======================" << std::endl;
	std::cout << "   System Parameters" << std::endl;
	std::cout << "=======================" << std::endl;
	std::cout << "n_particles =     " << n_particles << std::endl;
	std::cout << "number density =  " << number_density << std::endl;
	std::cout << "Lx =  " << L[0] << std::endl;
	std::cout << "Ly =  " << L[1] << std::endl;
	std::cout << "Lz =  " << L[2] << std::endl;
	std::cout << "=======================" << std::endl << std::endl;;
	
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
	std::ofstream dataOut("relax.xyz");	
    print_xyz(dataOut, x, L);

	// Number of accepted moves
	int n_accept=0; 

	// System potential
	double potential =0;	
	
	// Calculate the initial potential of the system
	calc_pe(x, nbr, sigma, epsilon, cutoff, L, &potential);
	std::cout << "Initial potential =  " << potential << std::endl;

	/*
        ===============================================
        Time loop to relax from initial configuration.
        Also, we will use this loop to find the optimal
        displacement for a given acceptance probability.
        ===============================================
    */
	
	for(int time = 0; time < n_relax; time++){
		
		// Number of moves accepted during each timestep
		int n_accept_local;
		
		// Attempt to translate each of the particles in the system, return
        // the percentage of moves accepted
		n_accept_local = pair_LJ_neighborlist(x, nbr, sigma, epsilon, cutoff, L,
                                              dx, &potential, beta);

		// Update the total number of accepted moves
		n_accept += n_accept_local; 
		
		// Check to see if we need to update our neighborlist
		nsq_neighbor_check(x, nbr, cutoff, skin, L);
		
		// Check to see if we should write data to file and to the screen
		if(time % adjust_freq == 0 && time > 0){

			// Calculate the total probability of accepted moves
			double prob = (double)n_accept / (double)(adjust_freq * n_particles);
			
			// Scale up or down the max displacement to achieve the
            // acceptance probability
			if(prob < target)
				dx = dx / 1.01;
			else if(prob > target)
				dx = dx * 1.01;
			
			if(dx > skin / 2.0)
				dx = skin / 2.0;

			std::cout << "Step: " << time << " of " << n_relax << "\tdx: " << dx << "\tprobability: " << prob << "\ttarget prob: " << target << std::endl;
			
			// Reset total acceptance accumulator 
			n_accept = 0;

			// Write trajectory data
			print_xyz(dataOut, x, L);
		}
	}
    std::cout << "Optimal displacement: " << dx << std::endl;
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
	double target = atof(argv[7]);
    int seed = atoi(argv[8]);
    int n_relax = atoi(argv[9]);
    int adjust_freq = atoi(argv[10]);

	relax(sigma, epsilon, cutoff, "system_init.gro", temperature, dx, target,
          seed, n_relax, adjust_freq);
	
	return 0;
}
