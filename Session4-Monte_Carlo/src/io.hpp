/*
 *  io.hpp
 *
 *  Created by Christopher Iacovella on 8/1/11.
 *  Modified by Andrew Z. Summers on 6/18/18.
 *  Copyright 2011 . All rights reserved.
 *
 */


std::ostream& operator << (std::ostream& os, coordlist_t& x)
{
	for (unsigned int i=0; i<x.size(); i++) {
		os << x[i] <<"\n";
	}
	return os;
}

std::ostream& operator << (std::ostream& os, coord_t& x)
{
	unsigned int sm1 = x.size()-1;
	for (unsigned int k=0; k<x.size(); k++) {
		os << x[k];
		if (k != sm1) {
			os <<"\t";
		}
	}
	return os;
}

std::istream& operator >> (std::istream& is, coordlist_t& x)
{
	return is;
}

std::istream& operator >> (std::istream& is, coord_t& x)
{
	return is;
}


void print(coordlist_t& x, std::ostream& os) 
{
	os << x << "\n";
}

void print(coord_t& x, std::ostream& os)
{
	os << x << "\n";
}

void load_raw(const char* filename, coordlist_t& x)
{
	std::ifstream file(filename);
	if (file.fail()) {
		std::cerr << "Error: load: can't open file " 
		<< filename << ".\n";
		exit(1);
	}
	
	std::string str;
	while (std::getline(file, str)) {
		std::istringstream iss(str.c_str());
		double temp;
		coord_t xi;
		while (iss >> temp) {
			xi.push_back(temp);
		}
		if (xi.size() != 0) {
			x.push_back(xi);
		}
	}
} 

void print_xyz(const char* filename, coordlist_t& x, double *L)
{
	std::ofstream dataOut(filename);	
	
	print_xyz(dataOut,x,L);

}

void print_xyz(std::ostream& dataOut, coordlist_t& x, double *L)
{	
	dataOut << x.size() << std::endl;
	dataOut << "#\t" << L[0] << "\t" << L[1] << "\t" << L[2] << std::endl;
	for(int i=0; i<x.size(); i++)
	{
		dataOut << "N\t" << x[i] << std::endl;
	}
	
}

void load_gro(const char* filename, coordlist_t& xyz, double *L)
{
    std::ifstream lp(filename);
    double x, y, z;
    int n_particles;
    std::string particle_type;
    
    // Make sure we can open the file
    ASSERT(lp.good(), "File " << filename << " not found.");
    
    xyz.clear();
    std::string temp_string, param;
    getline(lp, temp_string);
    getline(lp, temp_string);
    std::istringstream iss(temp_string.c_str());
    iss >> n_particles;

    ASSERT(n_particles > 0, "System contains no particles!");
    
    int temp_index;
    std::string temp_element;
    std::string temp_residue;
    
    for(int i = 0; i < n_particles; i++)
    {
        getline(lp, temp_string);
        std::istringstream iss2(temp_string.c_str());
        
        iss2 >> temp_residue >> temp_element >> temp_index >> x >> y >> z;
        //std::cout << "RES: " << temp_residue << "\tELE: " << temp_element << "\tID: " << temp_index << "\tX: " << x << "\tY: " << y << "\tZ: " << z << std::endl;
        
        coord_t ctemp;
        ctemp.push_back(x);
        ctemp.push_back(y);
        ctemp.push_back(z);
        xyz.push_back(ctemp);
    }
    getline(lp, temp_string);
    std::istringstream iss3(temp_string.c_str());
    iss3 >> x >> y >> z;
    //std::cout << "BOXX: " << x << "\tBOXY: " << y << "\tBOXZ" << z << std::endl;
    L[0] = x;
    L[1] = y;
    L[2] = z;
}
