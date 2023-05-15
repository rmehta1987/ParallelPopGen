/*
 * run.cu
 *
 *      Author: David Lawrie
 */

#include "go_fish.cuh"
#include "spectrum.h"
#include <chrono>
#include <cstdlib>
#include <fstream>


auto run_mse_robertson_model(float num_sites, float sel_coef, float mu, unsigned int seed1, unsigned int seed2){


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
    cp(3372), 54188,
    cp(3287), 54527,
    cp(3359), 54688,
    cp(3570), 54844,
    cp(4713), 55142,
    cp(5661), 55284,
    cp(7540), 55423,
    cp(11375), 55557,
    cp(14310), 55688,
    cp(13292), 55816,
    cp(14522), 54361,
    cp(613285), 55890,
    cp(3000000), 55940); // last one is the additional 50 generations 
	return GO_Fish::run_sim(input, mutation_rate , demography_model, migration, Sim_Model::make_stabilizing_cselection_model(sel_coef, 0.0f), inbreeding, Sim_Model::make_stabilizing_cdominance_model(), Sampling::off(), GO_Fish::allele_trajectories(), Sim_Model::stabilizing_mse_integrand());
}
void print_mse_robertson_sfs(int num_iterations, int num_samples, float sel_coef, float num_sites, float mu, std::string file_name)
{
	auto num_iter = num_iterations;															//number of iterations
	auto sample_size = num_samples;													//number of samples in SFS
    Spectrum::SFS my_spectra_mse;
    //std::vector<std::vector<float>> results(num_iter,  std::vector<float> (sample_size, 0)); // storage for SFS results
    std::vector<float> average(sample_size, 0.0);
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    GO_Fish::allele_trajectories b;
    
    std::cout << "Starting iteration 0 of simulation" << std::endl;
	for(auto j = 0; j < num_iterations; j++)
  {
		b = run_mse_robertson_model(num_sites, sel_coef, mu, 0xbeeff00d + 2 * j+3, 0xdecafbad - 2 * j+6);
		Spectrum::site_frequency_spectrum(my_spectra_mse,b,0,0,sample_size);
        for (int i = 0; i < sample_size; i++)
        {
            //results[j][i] = my_spectra.frequency_spectrum[i];
            // calculates average, // [ old_average * (n-1) + new_value ] / n
            if (j == 0)
            {
                average[i] += my_spectra_mse.frequency_spectrum[i]; 
            }
            else
            {
                average[i] = (average[i] * (j-1.0) + my_spectra_mse.frequency_spectrum[i]) * 1.0/j;
            }
        }
        if (j%2==0)
        {
             std::cout << "Finished iteration " << j << " of simulation" << std::endl;
        }
	}
  std::cout << "Finished all iterations of simulation" << std::endl;
	
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double,std::milli> elapsed_ms = end - start;
	std::cout<<"\nnumber of mutations in simulation: " << b.maximal_num_mutations() << "\nnumber of mutations in SFS: "<< my_spectra_mse.num_mutations <<"\ntime elapsed (ms): "<< elapsed_ms.count()/num_iter << std::endl << std::endl;

  std::ofstream output_file(file_name);
    //std::string output_file_name_debug = "last_sim_" + file_name;
    //std::ofstream output_file2(output_file_name_debug);
    std::cout<<"\n Saving SFS to: "<< file_name << std::endl;
    for (int i = 0; i < sample_size; i++)
    {
        output_file << average[i] << "\n";
        
        //output_file2 << my_spectra.frequency_spectrum[i] << "\n";;

    }

}

int main(int argc, char **argv)

{
    // this is a point selection coefficient the selection coefficient will remain the same for the population, this is the un-scaled selection coefficient
    float PointSel = .005;
    int num_samples = 1000; // Number of samples in gnomad: 111710
    std::string file_name = "out_file.txt";
    float num_sites = 133.0f * pow(10.f, 6);
    float mu = 1.24f * pow(10.f, -8); // genome-wide mutation rate https://www.biorxiv.org/content/10.1101/2022.07.11.499645v1.full.pdf
    int num_iterations = 1;
    if (argc != 4) // 3 Total parameters, [executable, unscaled selection coefficient, num_samples, file_seed]
    {
        fprintf(stderr, "Warning: The number of arguments given in the command line is not correct. In this version you need to pass in a selection cofficient, sample_size, and output_file_name, format is: ./GOFish scaled_mutation_rate unscaled_selection coefficient num_samples output_filename\n");
        // exit(8);
        std::cout << "Using default values" << std::endl;
    }
    else
    {
        PointSel = atof(argv[1]);
        num_samples = atoi(argv[2]);
        file_name = argv[3];

    }

    std::cout << "Currently we are using the genome wide mutation 1.5e-8: " << std::endl;
    std::cout << "Unscaled Point Selection: " << PointSel << std::endl;
    std::cout << "Number of samples to generate SFS: " << num_samples << std::endl;

    std::cout << "Running simulations" << std::endl;
    print_mse_robertson_sfs(num_iterations, num_samples, PointSel, num_sites, mu, file_name);
}
