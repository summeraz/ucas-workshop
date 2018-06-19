/*
 *  neighbor.cpp
 *  simple_LJ
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

void nsq_neighbor_init(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L) 
{
	nbr.clear();
	neighbor temp;
	for(int i=0; i< x.size(); i++)
	{
		nbr.push_back(temp);
	}
	double total_range = (cutoff + skin)*(cutoff + skin);

	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
	{
		L_inv[k] = 1.0/L[k];
	}
	for(int i=0; i< x.size(); i++)
	{
		nbr[i].member.push_back(i);
		nbr[i].x_old.push_back(x[i]);
		
		for(int j=0; j< x.size(); j++)
		{
			if(i!=j)
			{
				double r2=0;
				double rij[3];
				for(int k=0; k<3; k++)
				{
					rij[k] = x[i][k] - x[j][k];
					double pbc  = L[k]*anint(rij[k]*L_inv[k]);
					rij[k] -= pbc;
					r2 += rij[k]*rij[k];
				}
				if(r2 < total_range)
				{
					nbr[i].member.push_back(j);
					nbr[i].x_old.push_back(x[j]);	
			}		
			
			}
		}
	}
}

void nsq_neighbor_check(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L) 
{
	double half_skin2 = skin*skin/4.0;
	bool rebuild = false;
	
	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
	{
		L_inv[k] = 1.0/L[k];
	}
	for(int i=0; i<nbr.size(); i++)
	{
		if(rebuild == false)
		{
			for(int j=0; j<nbr[i].x_old.size(); j++)
			{
				int nbr_id = nbr[i].member[j];
				double r2=0;
				double rij[3];				
				for(int k=0; k<3; k++)
				{
					rij[k] = nbr[i].x_old[j][k] - x[nbr_id][k];
					double pbc  = L[k]*anint(rij[k]*L_inv[k]);
					rij[k] -= pbc;
					r2 += rij[k]*rij[k];
				}
				
				if(r2 > half_skin2)
				{
					rebuild = true;
					break;
				}
			}
		}
	}
	if(rebuild == true)
	{
		//std::cout << "need to rebuild" << std::endl;
		nsq_neighbor_rebuild(x, nbr, cutoff, skin, L);
	}
		
}
void nsq_neighbor_rebuild(coordlist_t& x, std::vector<neighbor> &nbr, double cutoff, double skin, double *L) 
{
	
	double total_range = (cutoff + skin)*(cutoff + skin);
	
	//precalculate 1/L
	double L_inv[3];
	for(int k=0; k<3; k++)
		L_inv[k] = 1.0/L[k];
	
	for(int i=0; i< x.size(); i++)
	{
		nbr[i].member.clear();
		nbr[i].x_old.clear();
		
		nbr[i].member.push_back(i);
		nbr[i].x_old.push_back(x[i]);
		
		for(int j=0; j< x.size(); j++)
		{
			if(i != j)
			{
				double r2=0;
				double rij[3];
				for(int k=0; k<3; k++)
				{
					rij[k] = x[i][k] - x[j][k];
					double pbc  = L[k]*anint(rij[k]*L_inv[k]);
					rij[k] -= pbc;
					r2 += rij[k]*rij[k];
				}
				if(r2 < total_range)
				{
					nbr[i].member.push_back(j);
					nbr[i].x_old.push_back(x[j]);	
		
				}
			}
			
		}
	}
}

