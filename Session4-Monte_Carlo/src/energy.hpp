/*
 *  force.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 1/25/12.
 *  Copyright 2012. All rights reserved.
 *
 */


int pair_LJ_neighborlist(coordlist_t& x, std::vector<neighbor> &nbr, double sigma, double epsilon, double cutoff, double *L, double dx, double *potential_energy_total, double _beta)
{
	
	
	//This function calculates the LJ potential between particles pairs of the form:
	//LJ_potential = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 )
	//LJ_force = 24*epsilon*( 2*(sigma^12)/(r^13) - (sigma^6)/(r^7) )
	
	
	//pre-calculate sigma^6 and sigma^12 for efficiency
	double sigma_p2 = sigma*sigma;
	double sigma_p6 = sigma_p2*sigma_p2*sigma_p2;
	double sigma_p12 = sigma_p6*sigma_p6;
	
	//precalculate the epsilon prefactors
	double epsilon_t4 = 4.0*epsilon;
	
	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
	
	double cutoff2 = cutoff*cutoff;
	double potential_energy = 0.0;
	
	//do not allow any particles to get closer than 0.85 to avoid overflow
	double inner_cutoff2 = 0.85*0.85;
	
	int n_accept=0;
	for(int i=0; i<x.size(); i++)
	{
		
	
		
		//calculate the potential energy of particle "i" with all of its neighbors "j"
		double potential_energy_old = 0.0;
		for(int j=1; j<nbr[i].member.size(); j++)
		{
			int nbr_id = nbr[i].member[j];
			
			//calculate the separation between two particles, 
			//taking into account periodic boundary conditions
			double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] = x[i][k] - x[nbr_id][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
			

			if(r2 < cutoff2)
			{
				//double r2inv = 1.0/r2;
				double r = sqrt(r2);
				double r_inv = 1.0/r;
				
				//precalculate a bunch of terms for efficiency
				double r_inv_p2 = r_inv*r_inv;			// 1/r^2
				double r_inv_p3 = r_inv*r_inv_p2;		// 1/r^3
				double r_inv_p6 = r_inv_p3*r_inv_p3;	// 1/r^6
				double r_inv_p12 = r_inv_p6*r_inv_p6;	// 1/r^12
				
				
				double LJ_potential = epsilon_t4*(sigma_p12*r_inv_p12 - sigma_p6*r_inv_p6);
				potential_energy_old += LJ_potential;
			}
			
		}

		//generate a random displacement of a particle
		//then check to see if we should accept this move
		coord_t new_x;
		for( int k=0; k<3; k++)
		{
			double dx_temp = dx*(2.0*drand48()-1.0);
			new_x.push_back(x[i][k]+dx_temp);
		}
		
		
		
		bool overlap = false;
		double potential_energy_new = 0.0;
		
		//calculate the potential energy of the translated particle "i" and all neighbors "j"
		for(int j=1; j<nbr[i].member.size(); j++)
		{
			int nbr_id = nbr[i].member[j];
			
			//calculate the separation between two particles, 
			//taking into account periodic boundary conditions
			double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] = new_x[k] - x[nbr_id][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
			
			//we don't want any particles less than 0.85
			//if r < 0.85, we automatically reject this move
			if(r2 < inner_cutoff2)
			{
				overlap=true;
				break;
			}
			else if(r2 < cutoff2)
			{
				//double r2inv = 1.0/r2;
				double r = sqrt(r2);
				double r_inv = 1.0/r;
				
				//precalculate a bunch of terms for efficiency
				double r_inv_p2 = r_inv*r_inv;			// 1/r^2
				double r_inv_p3 = r_inv*r_inv_p2;		// 1/r^3
				double r_inv_p6 = r_inv_p3*r_inv_p3;	// 1/r^6
				double r_inv_p12 = r_inv_p6*r_inv_p6;	// 1/r^12
				
		
				double LJ_potential = epsilon_t4*(sigma_p12*r_inv_p12 - sigma_p6*r_inv_p6);
				potential_energy_new += LJ_potential;

			}
		}
		
		if(overlap == false)
		{
			//calculate the change in potential energy related to the displacement
			double delta_PE = potential_energy_new - potential_energy_old;
			
			//if we lower the potential energy always accept
			if(delta_PE < 0)
			{
				*potential_energy_total = *potential_energy_total + delta_PE; 
				n_accept++;
				
				//assign the translated position to the main position array
				//making sure to apply PBC
				for( int k=0; k<3; k++)
				{
					new_x[k] -= L[k]*anint(new_x[k]*L_inv[k]);
					x[i][k] = new_x[k];
					
				}
			}
			else
			{
				//if the potential energy increases, check to see if we should accept
				double rand_value = drand48();
				if(exp(-_beta*delta_PE) > rand_value)
				{
					*potential_energy_total = *potential_energy_total + delta_PE; 
					n_accept++;
					
					//assign the translated position to the main position array
					//making sure to apply PBC
					for( int k=0; k<3; k++)
					{
						new_x[k] -= L[k]*anint(new_x[k]*L_inv[k]);
						x[i][k] = new_x[k];
						
					}
					
				}
			}
		}
		
	}
	return n_accept;
	
	
}

void calc_pe(coordlist_t& x, std::vector<neighbor> &nbr, double sigma, double epsilon, double cutoff, double *L, double *potential_energy_total)
{
	
	
	//This function calculates the LJ potential between particles pairs of the form:
	//LJ_potential = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 )
	//LJ_force = 24*epsilon*( 2*(sigma^12)/(r^13) - (sigma^6)/(r^7) )
	
	
	//pre-calculate sigma^6 and sigma^12 for efficiency
	double sigma_p2 = sigma*sigma;
	double sigma_p6 = sigma_p2*sigma_p2*sigma_p2;
	double sigma_p12 = sigma_p6*sigma_p6;
	
	//precalculate the epsilon prefactors
	double epsilon_t4 = 4.0*epsilon;
	
	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
	
	double cutoff2 = cutoff*cutoff;
	
	double potential_energy = 0.0;
	
	//brute force potential energy calculation, as it is done infrequently and does not need to be optimized
	for(int i=0; i<x.size(); i++)
	{
						
		for(int j=i+1; j<x.size(); j++)
		{
			
			//calculate the separation between two particles, 
			//taking into account periodic boundary conditions
			double rij[3];
			double r2=0;
			for(int k=0; k<3; k++)
			{
				rij[k] =x[i][k] - x[j][k];
				double pbc  = L[k]*anint(rij[k]*L_inv[k]);
				rij[k] -= pbc;
				r2 += rij[k]*rij[k];
			}
			
			if(r2 < cutoff2)
			{
				//double r2inv = 1.0/r2;
				double r = sqrt(r2);
				double r_inv = 1.0/r;
				
				//precalculate a bunch of terms for efficiency
				double r_inv_p2 = r_inv*r_inv;			// 1/r^2
				double r_inv_p3 = r_inv*r_inv_p2;		// 1/r^3
				double r_inv_p6 = r_inv_p3*r_inv_p3;	// 1/r^6
				double r_inv_p12 = r_inv_p6*r_inv_p6;	// 1/r^12
				
				
				double LJ_potential = epsilon_t4*(sigma_p12*r_inv_p12 - sigma_p6*r_inv_p6);
				potential_energy += LJ_potential;
			}
			
		}
		
	}
	*potential_energy_total = potential_energy;
	
}
