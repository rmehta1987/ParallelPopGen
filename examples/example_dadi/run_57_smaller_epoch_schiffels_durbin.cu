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
    
    0 14448
    4545 14068
    11956 14464
    17873 15208
    22799 16256
    27018 17618
    30709 19347
    33989 21534
    36940 24236
    39622 27367
    42081 30416
    44350 32060
    46458 31284
    47457 29404
    48424 26686
    49360 23261
    50268 18990
    50420 16490
    51123 12958
    51737 9827
    52284 7477
    52775 5791
    53222 4670
    53632 3841
    54010 3372
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

      */

void run_validation_test(const float sel_coef, const int num_samples, const std::string file_name)
{
    typedef std::vector<Sim_Model::demography_constant> dem_const;                                                        // vector of all the population sizes
    typedef Sim_Model::demography_piecewise<Sim_Model::demography_constant, Sim_Model::demography_constant> epoch_0_to_1; // initiali expansion
    typedef Sim_Model::demography_piecewise<epoch_0_to_1, Sim_Model::demography_constant> epoch_1_to_2;
    typedef Sim_Model::demography_piecewise<epoch_1_to_2, Sim_Model::demography_constant> epoch_2_to_3;
    typedef Sim_Model::demography_piecewise<epoch_2_to_3, Sim_Model::demography_constant> epoch_3_to_4;
    typedef Sim_Model::demography_piecewise<epoch_3_to_4, Sim_Model::demography_constant> epoch_4_to_5;
    typedef Sim_Model::demography_piecewise<epoch_4_to_5, Sim_Model::demography_constant> epoch_5_to_6;
    typedef Sim_Model::demography_piecewise<epoch_5_to_6, Sim_Model::demography_constant> epoch_6_to_7;
    typedef Sim_Model::demography_piecewise<epoch_6_to_7, Sim_Model::demography_constant> epoch_7_to_8;
    typedef Sim_Model::demography_piecewise<epoch_7_to_8, Sim_Model::demography_constant> epoch_8_to_9;
    typedef Sim_Model::demography_piecewise<epoch_8_to_9, Sim_Model::demography_constant> epoch_9_to_10;
    typedef Sim_Model::demography_piecewise<epoch_9_to_10, Sim_Model::demography_constant> epoch_10_to_11;
    typedef Sim_Model::demography_piecewise<epoch_10_to_11, Sim_Model::demography_constant> epoch_11_to_12;
    typedef Sim_Model::demography_piecewise<epoch_11_to_12, Sim_Model::demography_constant> epoch_12_to_13;
    typedef Sim_Model::demography_piecewise<epoch_12_to_13, Sim_Model::demography_constant> epoch_13_to_14;
    typedef Sim_Model::demography_piecewise<epoch_13_to_14, Sim_Model::demography_constant> epoch_14_to_15;
    typedef Sim_Model::demography_piecewise<epoch_14_to_15, Sim_Model::demography_constant> epoch_15_to_16;
    typedef Sim_Model::demography_piecewise<epoch_15_to_16, Sim_Model::demography_constant> epoch_16_to_17;
    typedef Sim_Model::demography_piecewise<epoch_16_to_17, Sim_Model::demography_constant> epoch_17_to_18;
    typedef Sim_Model::demography_piecewise<epoch_17_to_18, Sim_Model::demography_constant> epoch_18_to_19;
    typedef Sim_Model::demography_piecewise<epoch_18_to_19, Sim_Model::demography_constant> epoch_19_to_20;
    typedef Sim_Model::demography_piecewise<epoch_19_to_20, Sim_Model::demography_constant> epoch_20_to_21;
    typedef Sim_Model::demography_piecewise<epoch_20_to_21, Sim_Model::demography_constant> epoch_21_to_22;
    typedef Sim_Model::demography_piecewise<epoch_21_to_22, Sim_Model::demography_constant> epoch_22_to_23;
    typedef Sim_Model::demography_piecewise<epoch_22_to_23, Sim_Model::demography_constant> epoch_23_to_24;
    typedef Sim_Model::demography_piecewise<epoch_23_to_24, Sim_Model::demography_constant> epoch_24_to_25;
    typedef Sim_Model::demography_piecewise<epoch_24_to_25, Sim_Model::demography_constant> epoch_25_to_26;
    typedef Sim_Model::demography_piecewise<epoch_25_to_26, Sim_Model::demography_constant> epoch_26_to_27;
    typedef Sim_Model::demography_piecewise<epoch_26_to_27, Sim_Model::demography_constant> epoch_27_to_28;
    typedef Sim_Model::demography_piecewise<epoch_27_to_28, Sim_Model::demography_constant> epoch_28_to_29;
    typedef Sim_Model::demography_piecewise<epoch_28_to_29, Sim_Model::demography_constant> epoch_29_to_30;
    typedef Sim_Model::demography_piecewise<epoch_29_to_30, Sim_Model::demography_constant> epoch_30_to_31;
    typedef Sim_Model::demography_piecewise<epoch_30_to_31, Sim_Model::demography_constant> epoch_31_to_32;
    typedef Sim_Model::demography_piecewise<epoch_31_to_32, Sim_Model::demography_constant> epoch_32_to_33;
    typedef Sim_Model::demography_piecewise<epoch_32_to_33, Sim_Model::demography_constant> epoch_33_to_34;
    typedef Sim_Model::demography_piecewise<epoch_33_to_34, Sim_Model::demography_constant> epoch_34_to_35;
    typedef Sim_Model::demography_piecewise<epoch_34_to_35, Sim_Model::demography_constant> epoch_35_to_36;
    typedef Sim_Model::demography_piecewise<epoch_35_to_36, Sim_Model::demography_constant> epoch_36_to_37;
   

    typedef Sim_Model::migration_constant_equal mig_const; // no migration

    GO_Fish::allele_trajectories b;
    b.sim_input_constants.num_populations = 1; // number of populations
    b.sim_input_constants.num_generations = 60000;
    b.sim_input_constants.num_sites = 36.0f * pow(10.f, 6); // Should be 36 Megabase pairs
    b.sim_input_constants.compact_interval = 5;
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
    pop_history.push_back(14448);
    pop_history.push_back(14068);
    pop_history.push_back(14464);
    pop_history.push_back(15208);
    pop_history.push_back(16256);
    pop_history.push_back(17618);
    pop_history.push_back(19347);
    pop_history.push_back(21534);
    pop_history.push_back(24236);
    pop_history.push_back(27367);
    pop_history.push_back(30416);
    pop_history.push_back(32060);
    pop_history.push_back(31284);
    pop_history.push_back(29404);
    pop_history.push_back(26686);
    pop_history.push_back(23261);
    pop_history.push_back(18990);
    pop_history.push_back(16490);
    pop_history.push_back(12958);
    pop_history.push_back(9827);
    pop_history.push_back(7477);
    pop_history.push_back(5791);
    pop_history.push_back(4670);
    pop_history.push_back(3841);
    pop_history.push_back(3372);
    pop_history.push_back(3287);
    pop_history.push_back(3359);
    pop_history.push_back(3570);
    pop_history.push_back(4095);
    pop_history.push_back(4713);
    pop_history.push_back(5661);
    pop_history.push_back(7540);
    pop_history.push_back(11375);
    pop_history.push_back(14310);
    pop_history.push_back(13292);
    pop_history.push_back(14522);
    pop_history.push_back(613285);
    pop_history.push_back(1000000); // final size is 1,000,000

    //for (int i = 0; i < pop_history.size(); i++)
    //{
    //    std::cout << "pop_history: " << pop_history[i].N << std::endl;
    //}

    // Time points from the when a population size change occurs
    std::cout << "Creating inflection points" << std::endl;
    inflection_points.push_back(0); // ignore this index usually just here for reference
    inflection_points.push_back(4545);
    inflection_points.push_back(11956);
    inflection_points.push_back(17873);
    inflection_points.push_back(22799);
    inflection_points.push_back(27018);
    inflection_points.push_back(30709);
    inflection_points.push_back(33989);
    inflection_points.push_back(36940);
    inflection_points.push_back(39622);
    inflection_points.push_back(42081);
    inflection_points.push_back(44350);
    inflection_points.push_back(46458);
    inflection_points.push_back(47457);
    inflection_points.push_back(48424);
    inflection_points.push_back(49360);
    inflection_points.push_back(50268);
    inflection_points.push_back(50420);
    inflection_points.push_back(51123);
    inflection_points.push_back(51737);
    inflection_points.push_back(52284);
    inflection_points.push_back(52775);
    inflection_points.push_back(53222);
    inflection_points.push_back(53632);
    inflection_points.push_back(54010);
    inflection_points.push_back(54361);
    inflection_points.push_back(54527);
    inflection_points.push_back(54688);
    inflection_points.push_back(54844);
    inflection_points.push_back(54995);
    inflection_points.push_back(55142);
    inflection_points.push_back(55284);
    inflection_points.push_back(55423);
    inflection_points.push_back(55557);
    inflection_points.push_back(55688);
    inflection_points.push_back(55816);
    inflection_points.push_back(55890);
    inflection_points.push_back(55940); // last 50 generations had a pop size of one million
    //for (int i = 0; i < pop_history.size(); i++)
    //{
    //    std::cout << "inflection_point: " << inflection_points[i] << std::endl;
    //}
    std::cout << "Pop_history and inflection point size, should be the same: " << pop_history.size() << " " << inflection_points.size() << std::endl;

    // Construct the different epochs
    std::cout << "Creating the demographic history of epochs" << std::endl;
    epoch_0_to_1 epoch_0(pop_history[0], pop_history[1], inflection_points[1]);
    epoch_1_to_2 epoch_1(epoch_0, pop_history[2], inflection_points[2]);
    epoch_2_to_3 epoch_2(epoch_1, pop_history[3], inflection_points[3]);
    epoch_3_to_4 epoch_3(epoch_2, pop_history[4], inflection_points[4]);
    epoch_4_to_5 epoch_4(epoch_3, pop_history[5], inflection_points[5]);
    epoch_5_to_6 epoch_5(epoch_4, pop_history[6], inflection_points[6]);
    epoch_6_to_7 epoch_6(epoch_5, pop_history[7], inflection_points[7]);
    epoch_7_to_8 epoch_7(epoch_6, pop_history[8], inflection_points[8]);
    epoch_8_to_9 epoch_8(epoch_7, pop_history[9], inflection_points[9]);
    epoch_9_to_10 epoch_9(epoch_8, pop_history[10], inflection_points[10]);
    epoch_10_to_11 epoch_10(epoch_9, pop_history[11], inflection_points[11]);
    epoch_11_to_12 epoch_11(epoch_10, pop_history[12], inflection_points[12]);
    epoch_12_to_13 epoch_12(epoch_11, pop_history[13], inflection_points[13]);
    epoch_13_to_14 epoch_13(epoch_12, pop_history[14], inflection_points[14]);
    epoch_14_to_15 epoch_14(epoch_13, pop_history[15], inflection_points[15]);
    epoch_15_to_16 epoch_15(epoch_14, pop_history[16], inflection_points[16]);
    epoch_16_to_17 epoch_16(epoch_15, pop_history[17], inflection_points[17]);
    epoch_17_to_18 epoch_17(epoch_16, pop_history[18], inflection_points[18]);
    epoch_18_to_19 epoch_18(epoch_17, pop_history[19], inflection_points[19]);
    epoch_19_to_20 epoch_19(epoch_18, pop_history[20], inflection_points[20]);
    epoch_20_to_21 epoch_20(epoch_19, pop_history[21], inflection_points[21]);
    epoch_21_to_22 epoch_21(epoch_20, pop_history[22], inflection_points[22]);
    epoch_22_to_23 epoch_22(epoch_21, pop_history[23], inflection_points[23]);
    epoch_23_to_24 epoch_23(epoch_22, pop_history[24], inflection_points[24]);
    epoch_24_to_25 epoch_24(epoch_23, pop_history[25], inflection_points[25]);
    epoch_25_to_26 epoch_25(epoch_24, pop_history[26], inflection_points[26]);
    epoch_26_to_27 epoch_26(epoch_25, pop_history[27], inflection_points[27]);
    epoch_27_to_28 epoch_27(epoch_26, pop_history[28], inflection_points[28]);
    epoch_28_to_29 epoch_28(epoch_27, pop_history[29], inflection_points[29]);
    epoch_29_to_30 epoch_29(epoch_28, pop_history[30], inflection_points[30]);
    epoch_30_to_31 epoch_30(epoch_29, pop_history[31], inflection_points[31]);
    epoch_31_to_32 epoch_31(epoch_30, pop_history[32], inflection_points[32]);
    epoch_32_to_33 epoch_32(epoch_31, pop_history[33], inflection_points[33]);
    epoch_33_to_34 epoch_33(epoch_32, pop_history[34], inflection_points[34]);
    epoch_34_to_35 epoch_34(epoch_33, pop_history[35], inflection_points[35]);
    epoch_35_to_36 epoch_35(epoch_34, pop_history[36], inflection_points[36]);
    std::cout << "Starting final demography epoch: " << pop_history[37].N << " inflection point: " << inflection_points[37] << std::endl; 

    epoch_36_to_37 epoch_36(epoch_35, pop_history[37], inflection_points[37]); // final population size of one million
    std::cout << "Finished creating demographic events" << std::endl;

    // Migration parameters, no--migration
    mig_const mig_model;

    // Selection parameters
    Sim_Model::selection_constant weak_del((float)sel_coef);

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
        GO_Fish::run_sim(b, mutation, epoch_36, mig_model, weak_del, outbred, codominant, Sim_Model::bool_off(), Sim_Model::bool_off());
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
        fprintf(stderr, "Error: The number of arguments given in the command line is not correct. In this version you need to pass in a selection cofficient and unscaled mutation rate, format is: ./GOFish scaled_mutation_rate unscaled_selection coefficient num_samples \n");
        // exit(8);
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

    run_validation_test(PointSel, num_samples, file_name);
}
