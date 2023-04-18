/*
 * run.cu
 *
 *      Author: Rahul Mehta
 */

#include "go_fish.cuh"
#include "spectrum.h"
#include <vector>
#include <fstream>
#include <iterator>

/*
This folllows the simple zig-zag demographic model of schiffels-dubin for a single population:
https://github.com/popsim-consortium/demes-python/blob/main/examples/zigzag.yaml


description: A single population model with epochs of exponential growth and decay.
doi:
  - https://doi.org/10.1038/ng.3015
time_units: generations
demes:
  - name: generic
    description: All epochs wrapped into the same population, so that epoch intervals
      do not overlap, and they tile the entire existence of the population (all time,
      in this case).

epochs:
   0 - end_time: 34133.31
      start_size: 7156
   1 - end_time: 8533.33
      end_size: 71560
   2 - end_time: 2133.33
      end_size: 7156
   3 - end_time: 533.33
      end_size: 71560
   4 - end_time: 133.33
      end_size: 7156
   5 - end_time: 33.333
      end_size: 71560
   6 - end_time: 0
      end_size: 71560

I'm pretty sure these are going backwards (coalescent) as in the simulation starts 34133 generations ago and ends 0 generations ago (modern day), for forward simulations you need to invert the order of the generations e.g.:

epochs:
   0 - end_time: 0
      end_size: 7156
   1 - end_time: 25,599.98 (= 34133.31 - 8533.33 + 0)
      end_size: 71560
   2 - end_time: 31,999.98 (= 8533.33 - 2133.33 + 25,599.98)
      end_size: 7156
   3 - end_time: 33,599.98
      end_size: 71560
   4 - end_time: 33,999.98
      end_size: 7156
   5 - end_time: 34,099.98
      end_size: 71560
   6 - end_time: 34133.31
      end_size: 71560
 */

void run_validation_test(float mut_rate, float sel_coef, int num_samples){

    typedef Sim_Model::migration_constant_equal mig_const; // no migration
	

	GO_Fish::allele_trajectories b;
	b.sim_input_constants.num_populations = 1; 							//number of populations

    b.sim_input_constants.num_generations = 34150;
    b.sim_input_constants.num_sites = 36.f*pow(10.f,6);	 // Should be 36 Megabase pairs 
    b.sim_input_constants.compact_interval = 10;
    // Mutation and dominance parameters TODO Change dominance paramater to that of stabalizing selection

	Sim_Model::F_mu_h_constant codominant(0.5f); 						//dominance (co-dominant)
	Sim_Model::F_mu_h_constant outbred(0.f); 							//inbreeding (outbred)

	//Sim_Model::F_mu_h_constant mutation((float) mut_rate / (b.num_sites())); 	//per-site mutation rate 10^-9
    Sim_Model::F_mu_h_constant mutation(9.61111f*pow(10.f,-10)); 				//per-site 
    //std::cout << "mutation rate: " << 9.61f*pow(10.f,-10) << std::endl;
    

    // Demographic model
   	int N_ind = 7156;					//initial number of individuals in population
	Sim_Model::demography_constant demography_model(N_ind);	
 
    // Migration parameters, no--migration
    mig_const mig_model;

    // Selection parameters
	Sim_Model::selection_constant weak_del((float) sel_coef);

    // SFS parameters
	int sample_size = num_samples;										//number of samples in SFS
	int num_iter = 50;													//number of iterations
    Spectrum::SFS my_spectra;

    cudaEvent_t start, stop;											//CUDA timing functions
    float elapsedTime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	float avg_num_mutations = 0;
	float avg_num_mutations_sim = 0;
	std::vector<std::vector<float> > results(num_iter); 				//storage for SFS results
	for(int j = 0; j < num_iter; j++){ results[j].reserve(sample_size); }

	for(int j = 0; j < num_iter; j++){
		if(j == num_iter/2){ cudaEventRecord(start, 0); } 				//use 2nd half of the simulations to time simulation runs + SFS creation

		b.sim_input_constants.seed1 = 0xbeeff00d + 2*j; 				//random number seeds
		b.sim_input_constants.seed2 = 0xdecafbad - 2*j;
		GO_Fish::run_sim(b, mutation, demography_model, mig_model, weak_del, outbred, codominant, Sim_Model::bool_off(), Sim_Model::bool_off());
		Spectrum::site_frequency_spectrum(my_spectra,b,0,0,sample_size);

		avg_num_mutations += ((float)my_spectra.num_mutations)/num_iter;
		avg_num_mutations_sim += b.maximal_num_mutations()/num_iter;
        
		for(int i = 0; i < sample_size; i++){ results[j][i] = my_spectra.frequency_spectrum[i]; }
	}

	elapsedTime = 0;
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	//output SFS simulation results
	std::cout<<"SFS :"<<std::endl<< "allele count\tavg# mutations\tstandard dev\tcoeff of variation (aka relative standard deviation)"<< std::endl;
	for(int i = 1; i < sample_size; i++){
		double avg = 0;
		double std = 0;
		float num_mutations;
		for(int j = 0; j < num_iter; j++){ num_mutations = b.num_sites() - results[j][0]; avg += results[j][i]/(num_iter*num_mutations); }
		for(int j = 0; j < num_iter; j++){ num_mutations = b.num_sites() - results[j][0]; std += 1.0/(num_iter-1)*pow(results[j][i]/num_mutations-avg,2); }
		std = sqrt(std);
		std::cout<<i<<"\t"<<avg<<"\t"<<std<<"\t"<<(std/avg)<<std::endl;
	}

    std::ofstream output_file("./out_sfs.txt");
    for(int i = 0; i < sample_size; i++)
    { 
        output_file << my_spectra.frequency_spectrum[i] << "\n"; 
    }
    

	std::cout<<"\nnumber of sites in simulation: "<< b.num_sites() <<"\ncompact interval: "<< b.last_run_constants().compact_interval;
	std::cout<<"\naverage number of mutations in simulation: "<<avg_num_mutations_sim<<"\naverage number of mutations in SFS: "<<avg_num_mutations<<"\ntime elapsed (ms): "<< 2*elapsedTime/num_iter<<std::endl;
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

int main(int argc, char **argv) 

{ 
     // this is the mutation rate scaled with respect to number of sites, mutation_rate*(number of sites)
    float mut_rate = 0.3426;    
    // this is a point selection coefficient the selection coefficient will remain the same for the population, this is the un-scaled selection coefficient
    float PointSel = -.005; 
    int num_samples = 100;    

    // Number of samples for to generate for the site-frequency spectrum (SFS

    // Eventually this will read in a demographic history file for easier command line use instead of having to re-compile for every new demography <- possible but will still require a compilation step as Functors (functions passed as templates) need to be known at compile time (a requirement of GPUs), I have not yet added the ability to do this to the library, technically there are other libraries that will allow this, but I haven't merged them with my API to make it easy. It's on my TODO list


    if (argc != 4) // 3 Total parameters, [executable, scaled mutation rate, unscaled selection coefficient, num_samples]
    {
        fprintf(stderr, "Error: The number of arguments given in the command line is not correct. In this version you need to pass in a selection cofficient and unscaled mutation rate, format is: ./GOFish scaled_mutation_rate unscaled_selection coefficient num_samples \n");
        //exit(8);
        std::cout << "Using default values" << std::endl;
    }
    else{

        mut_rate = atof(argv[1]);
        PointSel = atof(argv[2]);
        num_samples = atoi(argv[3]);
    }

    std::cout<<"Scaled Mutation Rate: " << mut_rate << std::endl;
    std::cout<<"Inscaled Point Selection: " << PointSel << std::endl;
    std::cout<<"Number of samples to generate SFS: " << num_samples << std::endl;



    std::cout<<"Running simulations" << std::endl;

    run_validation_test(mut_rate, PointSel, num_samples); 
    
    }
