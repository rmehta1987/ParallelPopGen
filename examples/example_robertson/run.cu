/*
 * run.cu
 *
 *      Author: David Lawrie
 */

#include "go_fish.cuh"
#include "spectrum.h"
#include <chrono>
#include <cstdlib>

auto run_mse_robertson_model(float sel_coef, float variance, float N, float num_sites, float mu, unsigned int seed1, unsigned int seed2){




	// Simple parameters
	using cp = Sim_Model::constant_parameter;
  unsigned int num_gens = 55950;	    // Last Generation is 55940, so do another burn-in generations of 10
  unsigned int compact_interval = 10;
	GO_Fish::sim_constants input;
    Sim_Model::constant_parameter migration(0.0f); 							//constant migration rate
    //Sim_Model::constant_parameter selection(sel_coef); 				//constant, neutral, selection coefficient
	//Sim_Model::constant_parameter dominance(0.5f); 							//// dominance (co-dominant)
    Sim_Model::constant_parameter mutation_rate(mu); 			    // mutation rate
    Sim_Model::constant_parameter inbreeding(0.0f); 			    // no inbreeding
	//Sim_Model::constant_parameter variance(sel_coef); 		
    input.num_generations = num_gens;										//1,000 generations in simulation
	  input.seed1 = seed1; 													//random number seeds
	  input.seed2 = seed2;	
    input.num_sites = 	num_sites;            //number of sites
	  input.compact_interval = compact_interval;									//compact interval (in general: decrease compact interval for larger number of sites)

    // Complex UK Biobank Demographic Histroy
    //  I believe it goes - generation pop size_1, generation time_1 then eneration pop size_2, generation time_2
    
    auto demography_model = Sim_Model::make_piecewise_evolution_model(cp(14448), 4545,  
    cp(14068), 11956, 
    cp(14464), 15063,    
    cp(15208), 20439, 
    cp(16256), 24984, 
    cp(17618), 28922, 
    cp(19347), 32395, 
    cp(21534), 35501,
    cp(24236), 38312,
    cp(27367), 40877,
    cp(30416), 43238,
    cp(32060), 45423,
    cp(31284), 47457,
    cp(29404), 48424,
    cp(26686), 49360,
    cp(23261), 50268,
    cp(18990), 50420,
    cp(16490), 50784,
    cp(12958), 51440,
    cp(7477), 52775,
    cp(5791), 53004,
    cp(4670), 53431,
    cp(3841), 53824,
    //cp(3372), 54188,
    //cp(3287), 54527,
    //cp(3359), 54688,
    //cp(3570), 54844,
    //cp(4713), 55142,
    //cp(5661), 55284,
    //cp(7540), 55423,
    //cp(11375), 55557,
    //cp(14310), 55688,
    //cp(13292), 55816,
    cp(14522), 54361,
    cp(613285), 55890,
    cp(5000000), 55940); // last one is the additional 50 generations 
	return GO_Fish::run_sim(input, mutation_rate , demography_model, migration, Sim_Model::make_stabilizing_selection_model(sel_coef,variance), inbreeding, Sim_Model::make_stabilizing_dominance_model(), Sampling::off(), GO_Fish::allele_trajectories(), Sim_Model::stabilizing_mse_integrand());
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
	std::cout<<"\n"<< variance << ",";
	for(int i = 1; i < sample_size; i++){ std::cout<< my_spectra_mse.frequency_spectrum[i]/my_spectra_mse.num_mutations << ","; }
}

auto run_mse_model(float S, float N, float num_sites, float mu, unsigned int seed1, unsigned int seed2){
	using cp = Sim_Model::constant_parameter;
	Sim_Model::effective_parameter eff(N,0);
	return GO_Fish::run_sim({seed1,seed2,0,num_sites,1}, cp(mu), cp(N), cp(0), cp(eff(S)), cp(0), cp(0.5), Sampling::off());
}

void print_mse_sfs(int sample_size, float selection, float N, float num_sites, float mu){												
    Spectrum::SFS my_spectra_mse;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    GO_Fish::allele_trajectories b = run_mse_model(selection, N, num_sites, mu, 0xbeeff00d, 0xdecafbad);
	Spectrum::site_frequency_spectrum(my_spectra_mse,b,0,0,sample_size);

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double,std::milli> elapsed_ms = end - start;
	std::cout<<"\n"<< selection << ",";
	for(int i = 1; i < sample_size; i++){ std::cout<< my_spectra_mse.frequency_spectrum[i]/my_spectra_mse.num_mutations << ","; }
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

int main(int argc, char **argv) { 
	std::vector<float> variance = {-.005, -.01};
	//srand(0xbeeff00d);
	//std::vector<float> variance = {50, 50, 50, 50, 50, 50, 50, 50, 50};
	/*for(const auto & vs: variance){
		print_mse_robertson_sfs(51, vs, 0.5, 1*pow(10.f,4), 100*pow(10.f,7), 10*pow(10.f,-9)); 
	}*/
	print_mse_robertson_sfs(51, -.005, 0.5, 1*pow(10.f,4), 100*pow(10.f,7), 10*pow(10.f,-9)); 
	/*
	std::vector<float> selection = {25, 0, -10, -25, -85, -200};
	for(const auto & s: selection){
		print_mse_sfs(51, s, 1*pow(10.f,4), 100*pow(10.f,7), 10*pow(10.f,-9)); 
	}*/
}
