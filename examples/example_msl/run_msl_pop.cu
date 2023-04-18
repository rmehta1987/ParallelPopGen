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
//#include <thrust/device_vector.h>
//#include <thrust/host_vector.h>
//#include <eigen3/Eigen/Dense>

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

15105 10852
760 25030
0 66000
 */

void run_validation_test(const float sel_coef, const int num_samples, const std::string file_name)
{
	typedef std::vector<Sim_Model::demography_constant> dem_const;                                                        // vector of all the population sizes
	typedef Sim_Model::demography_piecewise<Sim_Model::demography_constant, Sim_Model::demography_constant> epoch_0_to_1; // initial expansion
	typedef Sim_Model::demography_exponential_growth exp_growth;
	typedef Sim_Model::demography_piecewise<epoch_0_to_1, exp_growth> demographic_model;
	typedef Sim_Model::demography_piecewise<epoch_0_to_1, Sim_Model::demography_constant> final_model; 
    
	//typedef Sim_Model::demography_piecewise<demographic_model, Sim_Model::demography_constant> final_model;


	typedef Sim_Model::migration_constant_equal mig_const; // no migration

	GO_Fish::allele_trajectories b;
	b.sim_input_constants.num_populations = 1; // number of populations
	b.sim_input_constants.num_generations = 15200;
	b.sim_input_constants.num_sites = 36.0f * pow(10.f, 6); // Should be 36 Megabase pairs
	b.sim_input_constants.compact_interval = 30;
	// Mutation and dominance parameters TODO Change dominance paramater to that of stabalizing selection

	Sim_Model::F_mu_h_constant codominant(0.5f); // dominance (co-dominant)
	Sim_Model::F_mu_h_constant outbred(0.f);     // inbreeding (outbred)

	// Sim_Model::F_mu_h_constant mutation((float) mut_rate / (b.num_sites())); 	//per-site mutation rate 10^-9
	Sim_Model::F_mu_h_constant mutation(9.6111f * pow(10.f, -9)); // per-site mutation rate -- testing -- need to get from command line
	//Sim_Model::F_mu_h_constant mutation(7.1111f * pow(10.f, -10)); // per-site mutation rate for the loss of function mutations

	// Demographic model
	// needs to be inverted as specified above
	std::vector<float> inflection_points;
	std::cout << "Creating population histories " << std::endl;
	// The different population sizes
	dem_const pop_history;
	pop_history.push_back(10852);
	pop_history.push_back(25030);
	pop_history.push_back(66000);

	std::cout << "Creating inflection points" << std::endl;
	inflection_points.push_back(0); // ignore this index usually just here for reference
	inflection_points.push_back(14345);
	inflection_points.push_back(15105);
	std::cout << "Pop_history and inflection point size, should be the same: " << pop_history.size() << " " << inflection_points.size() << std::endl;

	std::cout << "Creating the demographic history of epochs" << std::endl;
	
    epoch_0_to_1 epoch_0(pop_history[0], pop_history[1], inflection_points[1]);
	//exp_growth pop_exp_growth(0.12768f, 25030, 1);
	//demographic_model deme_final(epoch_0, pop_exp_growth, inflection_points[1]);
	final_model deme_final(epoch_0, pop_history[2], inflection_points[2]);

	std::cout << "Finished creating demographic events" << std::endl;

	// Migration parameters, no--migration
	mig_const mig_model;

	// Selection parameters
	//float gamma = -4;
	//Sim_Model::selection_constant weak_del(gamma, deme_final, outbred);
    Sim_Model::selection_constant weak_del(sel_coef);

	// SFS parameters
	int sample_size = num_samples; // number of samples in SFS
	const int num_iter = 5;              // number of iterations
	Spectrum::SFS my_spectra;

	cudaEvent_t start, stop; // CUDA timing functions
	float elapsedTime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	float avg_num_mutations = 0;
	float avg_num_mutations_sim = 0;
	std::vector<std::vector<float>> results(num_iter,  std::vector<float> (sample_size, 0)); // storage for SFS results
	std::vector<float> average(sample_size, 0.0);

	//thrust::device_vector<float> results; // ideally something like this eventually: https://stackoverflow.com/questions/21428378/reduce-matrix-columns-with-cuda
	std::cout << "Starting iteration 1 of simulation" << std::endl;
	for (int j = 0; j < num_iter; j++)
	{
		if (j == num_iter / 2)
		{
			cudaEventRecord(start, 0);
		} // use 2nd half of the simulations to time simulation runs + SFS creation

		b.sim_input_constants.seed1 = 0xbeeff00d + 2 * j; // random number seeds
		b.sim_input_constants.seed2 = 0xdecafbad - 2 * j;
		GO_Fish::run_sim(b, mutation, deme_final, mig_model, weak_del, outbred, codominant, Sim_Model::bool_off(), Sim_Model::bool_off());
		Spectrum::site_frequency_spectrum(my_spectra, b, 0, 0, sample_size);

		avg_num_mutations += ((float)my_spectra.num_mutations) / num_iter;
		avg_num_mutations_sim += b.maximal_num_mutations() / num_iter;

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
	}
	std::cout << "Finished all iterations of simulation" << std::endl;

	elapsedTime = 0;
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	std::ofstream output_file(file_name);
	//std::string output_file_name_debug = "last_sim_" + file_name;
	//std::ofstream output_file2(output_file_name_debug);

	for (int i = 0; i < sample_size; i++)
	{
		output_file << average[i] << "\n";
		//output_file2 << my_spectra.frequency_spectrum[i] << "\n";;

	}
/*
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

*/
	std::cout << "\nnumber of sites in simulation: " << b.num_sites() << "\ncompact interval: " << b.last_run_constants().compact_interval;
	std::cout << "\naverage number of mutations in simulation: " << avg_num_mutations_sim << "\naverage number of mutations in SFS: " << avg_num_mutations << "\ntime elapsed (ms): " << 2 * elapsedTime / num_iter << std::endl;
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

int main(int argc, char **argv)

{
	// this is a point selection coefficient the selection coefficient will remain the same for the population, this is the un-scaled selection coefficient
	float PointSel = -.005;
	int num_samples = 1000;
	std::string file_name = "out_file.txt";

	if (argc != 4) // 3 Total parameters, [executable, unscaled selection coefficient, num_samples, file_seed]
	{
		fprintf(stderr, "Error: The number of arguments given in the command line is not correct. In this version you need to pass in a selection cofficient, sample size, output file name, format is: ./GOFish unscaled_selection num_samples  file_name\n");
		// exit(8);
		std::cout << "Using default values" << std::endl;
	}
	else
	{
		PointSel = atof(argv[1]);
		num_samples = atoi(argv[2]);
		file_name = argv[3];

	}

	std::cout << "Currently we are using a scaled Mutation Rate pf .0256: " << std::endl;
	std::cout << "Unscaled Point Selection: " << PointSel << std::endl;
	std::cout << "Number of samples to generate SFS: " << num_samples << std::endl;

	std::cout << "Running simulations" << std::endl;

	run_validation_test(PointSel, num_samples, file_name);
}
