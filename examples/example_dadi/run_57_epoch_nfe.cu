/*
 * run_57_epoch_nfe.cu
 *
 *      Author: Rahul Mehta
 */

#include "../../3P/go_fish.cuh"
#include "spectrum.h"
#include <chrono>
#include <fstream>


/*
This folllows schiffels-durbin demographic model

time_units: generations
demes:
  - name: generic
    description: All epochs wrapped into the same population, so that epoch intervals
      do not overlap, and they tile the entire existence of the population (all time,
      in this case).

    0 14448
    4545 14068
    8483 14068
    11956 14464
    15063 14464
    17873 15208
    20439 15208
    22799 16256
    24984 16256
    27018 17618
    28922 17618
    30709 19347
    32395 19347
    33989 21534
    35501 21534
    36940 24236
    38312 24236
    39622 27367
    40877 27367
    42081 30416
    43238 30416
    44350 32060
    45423 32060
    46458 31284
    47457 29404
    48424 26686
    49360 23261
    50268 18990
    50420 16490
    50784 16490
    51123 12958
    51440 12958
    51737 9827
    52018 9827
    52284 7477
    52536 7477
    52775 5791
    53004 5791
    53222 4670
    53431 4670
    53632 3841
    53824 3841
    54010 3372
    54188 3372
    54361 3287
    54527 3359
    54688 3570
    54844 4095
    54995 4713
    55142 5661
    55284 7540
    55423 11375
    55557 14310
    55688 13292
    55816 14522
    55890 613285
    55940 5000000 -- for the last population growth over 50 years
*/

auto run_model(float num_sites, float mu, float sel_coef, unsigned int seed1, unsigned int seed2){

	using cp = Sim_Model::constant_parameter;
  unsigned int num_gens = 55950;	    // Last Generation is 55940, so do another burn-in generations of 10
  unsigned int compact_interval = 10;
  
  // Simple parameters
	GO_Fish::sim_constants input;
    Sim_Model::constant_parameter migration(0.0f); 							//constant migration rate
    Sim_Model::constant_parameter selection(sel_coef); 				//constant, neutral, selection coefficient
	Sim_Model::constant_parameter dominance(0.5f); 							//// dominance (co-dominant)
    Sim_Model::constant_parameter mutation_rate(mu); 			    // mutation rate
    Sim_Model::constant_parameter inbreeding(0.0f); 			    // no inbreeding
    input.num_generations = num_gens;										//1,000 generations in simulation
	  input.seed1 = seed1; 													//random number seeds
	  input.seed2 = seed2;	
    input.num_sites = 	num_sites;            //number of sites
	  input.compact_interval = compact_interval;									//compact interval (in general: decrease compact interval for larger number of sites)

    // Complex UK Biobank Demographic Histroy
    //  I believe it goes - generation pop size_1, generation time_1 then eneration pop size_2, generation time_2
    auto demography_model = Sim_Model::make_piecewise_evolution_model(cp(14448), 4545,  cp(14068), 8483, 
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
    cp(1000000), 55940); // last one is the additional 50 generations 
    

	return GO_Fish::run_sim(input, mutation_rate , demography_model, migration , selection, inbreeding, dominance, Sampling::off());
}

void print_sfs(int num_iterations, int num_samples, float num_sites, float mu, float sel_coef, std::string file_name){
	
  auto num_iter = num_iterations;															//number of iterations
	auto sample_size = num_samples;													//number of samples in SFS
  


    Spectrum::SFS my_spectra;
    GO_Fish::allele_trajectories b;
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<float>> results(num_iter,  std::vector<float> (sample_size, 0)); // storage for SFS results
    std::vector<float> average(sample_size, 0.0);

  std::cout << "Starting iteration 0 of simulation" << std::endl;
	for(auto j = 0; j < num_iter; j++)
  {
		b = run_model(num_sites, mu, sel_coef, 0xbeeff00d + 2 * j+3, 0xdecafbad - 2 * j+6);
		Spectrum::site_frequency_spectrum(my_spectra,b,0,0,sample_size);
    for (int i = 0; i < sample_size; i++)
        {
            results[j][i] = my_spectra.frequency_spectrum[i];
            // calculates average, // [ old_average * (n-1) + new_value ] / n
            if (j == 0)
            {
                average[i] += my_spectra.frequency_spectrum[i]; 
            }
            else
            {
                average[i] = (average[i] * (j-1.0) + my_spectra.frequency_spectrum[i]) * 1.0/j;
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
	std::cout<<"\nnumber of mutations in simulation: " << b.maximal_num_mutations() << "\nnumber of mutations in SFS: "<< my_spectra.num_mutations <<"\ntime elapsed (ms): "<< elapsed_ms.count()/num_iter << std::endl << std::endl;

  std::ofstream output_file(file_name);
    //std::string output_file_name_debug = "last_sim_" + file_name;
    //std::ofstream output_file2(output_file_name_debug);

    for (int i = 0; i < sample_size; i++)
    {
        output_file << average[i] << "\n";
        //output_file2 << my_spectra.frequency_spectrum[i] << "\n";;

    }
	
  /*for(int i = 1; i < sample_size; i++)
  { 
      std::cout<< my_spectra.frequency_spectrum[i]/my_spectra.num_mutations << std::endl; 
  }*/

  
    
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

int main(int argc, char **argv) 
{
    float num_sites = 36.0f * pow(10.f, 6); // Should be 36 Megabase pairs for gnomAD data, the largest chromosome (1) is 249 million nucleotide base pairs

    // this is a point selection coefficient the selection coefficient will remain the same for the population, this is the un-scaled selection coefficient
    float PointSel = -.005;
    int num_samples = 111710;
    float mutrate = (9.6111f * pow(10.f, -9));
    std::string file_name = "out_file.txt";

    // Number of samples for to generate for the site-frequency spectrum (SFS

    // Eventually this will read in a demographic history file for easier command line use instead of having to re-compile for every new demography <- possible but will still require a compilation step as Functors (functions passed as templates) need to be known at compile time (a requirement of GPUs), I have not yet added the ability to do this to the library, technically there are other libraries that will allow this, but I haven't merged them with my API to make it easy. It's on my TODO list

    if (argc != 3) // 3 Total parameters, [executable, unscaled selection coefficient, num_samples, file_name]
    {
        fprintf(stderr, "Warning: The number of arguments given in the command line is not correct. In this version you need to pass in a selection cofficient sample size, and output_file name, format is: ./GOFish unscaled_selection coefficient num_samples out_file_name\n");
        std::cout << "Using default values" << std::endl;
    }
    else
    {

        PointSel = atof(argv[1]);
        num_samples = atoi(argv[2]);
        file_name = argv[3];
    }

    std::cout << "Currently we are using a scaled Mutation Rate pf .3426 (missense mutation rate): " << std::endl;
    std::cout << "Unscaled Point Selection: " << PointSel << std::endl;
    std::cout << "Number of samples to generate SFS: " << num_samples << std::endl;

    std::cout << "Running simulations" << std::endl;

    print_sfs(2, num_samples, num_sites, mutrate, PointSel, file_name); 
     
    
}
