/*
 * run.cu
 *
 *      Author: David Lawrie
 */

#include "go_fish.cuh"
#include "spectrum.h"
#include <chrono>
#include <cstdlib>

auto run_mse_robertson_model(float effect_size, float variance, float N, float num_sites, float mu, unsigned int seed1, unsigned int seed2){
	using cp = Sim_Model::constant_parameter;
	return GO_Fish::run_sim({seed1,seed2,0,num_sites,1}, cp(mu), cp(N), cp(0), Sim_Model::make_stabilizing_cselection_model(effect_size,variance), cp(0), Sim_Model::make_stabilizing_cdominance_model(), Sampling::off(), GO_Fish::allele_trajectories(), Sim_Model::stabilizing_mse_integrand());
}

auto run_mse_robertson_model2(float effect_size, float variance, float N, float num_sites, float mu, unsigned int seed1, unsigned int seed2){
	using cp = Sim_Model::constant_parameter;
    return GO_Fish::run_sim({seed1,seed2,0,num_sites,1}, cp(mu), cp(N), cp(0), Sim_Model::make_stabilizing_selection_model(effect_size,variance), cp(0), Sim_Model::make_stabilizing_dominance_model(), Sampling::off(), GO_Fish::allele_trajectories(), Sim_Model::stabilizing_mse_integrand());
}


void print_mse_robertson_sfs(int sample_size, float effect_size, float variance, float N, float num_sites, float mu){												
    Spectrum::SFS my_spectra_mse;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    //GO_Fish::allele_trajectories b = run_mse_robertson_model(effect_size, variance, N, num_sites, mu, rand(), rand());
    GO_Fish::allele_trajectories b = run_mse_robertson_model(effect_size, variance, N, num_sites, mu, 0xbeeff00d, 0xdecafbad);
	Spectrum::site_frequency_spectrum(my_spectra_mse,b,0,0,sample_size);

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double,std::milli> elapsed_ms = end - start;
	std::cout<<"\n The effect sizeS is: "<< effect_size << "," << std::endl;
	for(int i = 1; i < sample_size; i++){ std::cout<< my_spectra_mse.frequency_spectrum[i]/my_spectra_mse.num_mutations << ","; }
    std::cout<<"\n Finished experiment with: "<< effect_size << std::endl;;
}

void print_mse_robertson_sfs2(int sample_size, float effect_size, float variance, float N, float num_sites, float mu){												
    Spectrum::SFS my_spectra_mse;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    //GO_Fish::allele_trajectories b = run_mse_robertson_model(effect_size, variance, N, num_sites, mu, rand(), rand());
    GO_Fish::allele_trajectories b = run_mse_robertson_model2(effect_size, variance, N, num_sites, mu, 0xbeeff00d, 0xdecafbad);
	Spectrum::site_frequency_spectrum(my_spectra_mse,b,0,0,sample_size);

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double,std::milli> elapsed_ms = end - start;
	std::cout<<"\n The variance is: "<< variance << ",";
	for(int i = 1; i < sample_size; i++){ std::cout<< my_spectra_mse.frequency_spectrum[i]/my_spectra_mse.num_mutations << ","; }
    std::cout<<"\n Finished experiment with: "<< variance << std::endl;;
}

auto run_mse_model(float S, float N, float num_sites, float mu, unsigned int seed1, unsigned int seed2){
	using cp = Sim_Model::constant_parameter;
	Sim_Model::effective_parameter eff(N,0);
	return GO_Fish::run_sim({seed1,seed2,0,num_sites,1}, cp(mu), cp(N), cp(0), cp(S), cp(0), cp(0.5), Sampling::off());
}

void print_mse_sfs(int sample_size, float selection, float N, float num_sites, float mu){												
    Spectrum::SFS my_spectra_mse;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    GO_Fish::allele_trajectories b = run_mse_model(selection, N, num_sites, mu, 0xbeeff00d, 0xdecafbad);
	Spectrum::site_frequency_spectrum(my_spectra_mse,b,0,0,sample_size);

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double,std::milli> elapsed_ms = end - start;
	std::cout<<"\n"<< selection << "," << std::endl;
	for(int i = 1; i < sample_size; i++){ std::cout<< my_spectra_mse.frequency_spectrum[i]/my_spectra_mse.num_mutations << ","; }
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

int main(int argc, char **argv) { 
    std::vector<float> es = {.04, 0.017, 0.01, 0.005, 0.002};
	//std::vector<float> variance = {0.05, 0, 0.03, 0.09, 0.45};
	//srand(0xbeeff00d);
	//std::vector<float> variance = {50, 50, 50, 50, 50, 50, 50, 50, 50};

    for(const auto & vs: es){
		print_mse_robertson_sfs(51, vs, 0.5, 1*pow(10.f,4), 1*pow(10.f,6), 10*pow(10.f,-8)); 
	}

    
    std::cout << "Doing original impelementation" << std::endl;
    std::vector<float> selection = {0, -.017, -.005, -.002};
	for(const auto & s: selection){
		print_mse_sfs(51, s, 1*pow(10.f,4), 1*pow(10.f,6), 10*pow(10.f,-8)); 
	}
    
}