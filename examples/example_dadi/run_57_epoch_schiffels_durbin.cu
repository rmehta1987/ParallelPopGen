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
*/

void run_validation_test(float mut_rate, float sel_coef, int num_samples)
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
    typedef Sim_Model::demography_piecewise<epoch_36_to_37, Sim_Model::demography_constant> epoch_37_to_38;
    typedef Sim_Model::demography_piecewise<epoch_37_to_38, Sim_Model::demography_constant> epoch_38_to_39;
    typedef Sim_Model::demography_piecewise<epoch_38_to_39, Sim_Model::demography_constant> epoch_39_to_40;
    typedef Sim_Model::demography_piecewise<epoch_39_to_40, Sim_Model::demography_constant> epoch_40_to_41;
    typedef Sim_Model::demography_piecewise<epoch_40_to_41, Sim_Model::demography_constant> epoch_41_to_42;
    typedef Sim_Model::demography_piecewise<epoch_41_to_42, Sim_Model::demography_constant> epoch_42_to_43;
    typedef Sim_Model::demography_piecewise<epoch_42_to_43, Sim_Model::demography_constant> epoch_43_to_44;
    typedef Sim_Model::demography_piecewise<epoch_43_to_44, Sim_Model::demography_constant> epoch_44_to_45;
    typedef Sim_Model::demography_piecewise<epoch_44_to_45, Sim_Model::demography_constant> epoch_45_to_46;
    typedef Sim_Model::demography_piecewise<epoch_45_to_46, Sim_Model::demography_constant> epoch_46_to_47;
    typedef Sim_Model::demography_piecewise<epoch_46_to_47, Sim_Model::demography_constant> epoch_47_to_48;
    typedef Sim_Model::demography_piecewise<epoch_47_to_48, Sim_Model::demography_constant> epoch_48_to_49;
    typedef Sim_Model::demography_piecewise<epoch_48_to_49, Sim_Model::demography_constant> epoch_49_to_50;
    typedef Sim_Model::demography_piecewise<epoch_49_to_50, Sim_Model::demography_constant> epoch_50_to_51;
    typedef Sim_Model::demography_piecewise<epoch_50_to_51, Sim_Model::demography_constant> epoch_51_to_52;
    typedef Sim_Model::demography_piecewise<epoch_51_to_52, Sim_Model::demography_constant> epoch_52_to_53;
    typedef Sim_Model::demography_piecewise<epoch_52_to_53, Sim_Model::demography_constant> epoch_53_to_54;
    typedef Sim_Model::demography_piecewise<epoch_53_to_54, Sim_Model::demography_constant> epoch_54_to_55;
    typedef Sim_Model::demography_piecewise<epoch_54_to_55, Sim_Model::demography_constant> epoch_55_to_56;

    typedef Sim_Model::migration_constant_equal mig_const; // no migration

    GO_Fish::allele_trajectories b;
    b.sim_input_constants.num_populations = 1; // number of populations

    b.sim_input_constants.num_generations = 34150;
    b.sim_input_constants.num_sites = 36.0f * pow(10.f, 6); // Should be 36 Megabase pairs
    b.sim_input_constants.compact_interval = 5;
    // Mutation and dominance parameters TODO Change dominance paramater to that of stabalizing selection

    Sim_Model::F_mu_h_constant codominant(0.5f); // dominance (co-dominant)
    Sim_Model::F_mu_h_constant outbred(0.f);     // inbreeding (outbred)

    // Sim_Model::F_mu_h_constant mutation((float) mut_rate / (b.num_sites())); 	//per-site mutation rate 10^-9
    Sim_Model::F_mu_h_constant mutation(9.6111f * pow(10.f, -9)); // per-site mutation rate -- testing -- need to get from command line 

    // Demographic model
    // needs to be inverted as specified above
    std::vector<float> inflection_points;
    std::cout << "Creating population histories " << std::endl;
    // The different population sizes
    dem_const pop_history;
    pop_history.push_back(14448);
    pop_history.push_back(14068);
    pop_history.push_back(14068);
    pop_history.push_back(14464);
    pop_history.push_back(14464);
    pop_history.push_back(15208);
    pop_history.push_back(15208);
    pop_history.push_back(16256);
    pop_history.push_back(16256);
    pop_history.push_back(17618);
    pop_history.push_back(17618);
    pop_history.push_back(19347);
    pop_history.push_back(19347);
    pop_history.push_back(21534);
    pop_history.push_back(21534);
    pop_history.push_back(24236);
    pop_history.push_back(24236);
    pop_history.push_back(27367);
    pop_history.push_back(27367);
    pop_history.push_back(30416);
    pop_history.push_back(30416);
    pop_history.push_back(32060);
    pop_history.push_back(32060);
    pop_history.push_back(31284);
    pop_history.push_back(29404);
    pop_history.push_back(26686);
    pop_history.push_back(23261);
    pop_history.push_back(18990);
    pop_history.push_back(16490);
    pop_history.push_back(16490);
    pop_history.push_back(12958);
    pop_history.push_back(12958);
    pop_history.push_back(9827);
    pop_history.push_back(9827);
    pop_history.push_back(7477);
    pop_history.push_back(7477);
    pop_history.push_back(5791);
    pop_history.push_back(5791);
    pop_history.push_back(4670);
    pop_history.push_back(4670);
    pop_history.push_back(3841);
    pop_history.push_back(3841);
    pop_history.push_back(3372);
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

    // Time points from the when a population size change occurs
    std::cout << "Creating inflection points" << std::endl;
    inflection_points.push_back(0); // ignore this index usually just here for reference
    inflection_points.push_back(4545);
    inflection_points.push_back(8483);
    inflection_points.push_back(11956);
    inflection_points.push_back(15063);
    inflection_points.push_back(17873);
    inflection_points.push_back(20439);
    inflection_points.push_back(22799);
    inflection_points.push_back(24984);
    inflection_points.push_back(27018);
    inflection_points.push_back(28922);
    inflection_points.push_back(30709);
    inflection_points.push_back(32395);
    inflection_points.push_back(33989);
    inflection_points.push_back(35501);
    inflection_points.push_back(36940);
    inflection_points.push_back(38312);
    inflection_points.push_back(39622);
    inflection_points.push_back(40877);
    inflection_points.push_back(42081);
    inflection_points.push_back(43238);
    inflection_points.push_back(44350);
    inflection_points.push_back(45423);
    inflection_points.push_back(46458);
    inflection_points.push_back(47457);
    inflection_points.push_back(48424);
    inflection_points.push_back(49360);
    inflection_points.push_back(50268);
    inflection_points.push_back(50420);
    inflection_points.push_back(50784);
    inflection_points.push_back(51123);
    inflection_points.push_back(51440);
    inflection_points.push_back(51737);
    inflection_points.push_back(52018);
    inflection_points.push_back(52284);
    inflection_points.push_back(52536);
    inflection_points.push_back(52775);
    inflection_points.push_back(53004);
    inflection_points.push_back(53222);
    inflection_points.push_back(53431);
    inflection_points.push_back(53632);
    inflection_points.push_back(53824);
    inflection_points.push_back(54010);
    inflection_points.push_back(54188);
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
    std::cout << "Finished 10 demographic events" << std::endl;
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
    std::cout << "Finished 20 demographic events" << std::endl;
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
    std::cout << "Finished 30 demographic events" << std::endl;
    epoch_30_to_31 epoch_30(epoch_29, pop_history[31], inflection_points[31]);
    epoch_31_to_32 epoch_31(epoch_30, pop_history[32], inflection_points[32]);
    epoch_32_to_33 epoch_32(epoch_31, pop_history[33], inflection_points[33]);
    epoch_33_to_34 epoch_33(epoch_32, pop_history[34], inflection_points[34]);
    epoch_34_to_35 epoch_34(epoch_33, pop_history[35], inflection_points[35]);
    epoch_35_to_36 epoch_35(epoch_34, pop_history[36], inflection_points[36]);
    epoch_36_to_37 epoch_36(epoch_35, pop_history[37], inflection_points[37]);
    epoch_37_to_38 epoch_37(epoch_36, pop_history[38], inflection_points[38]);
    epoch_38_to_39 epoch_38(epoch_37, pop_history[39], inflection_points[39]);
    epoch_39_to_40 epoch_39(epoch_38, pop_history[40], inflection_points[40]);
    std::cout << "Finished 40 demographic events" << std::endl;
        std::cout << "Population size at epoch 40: " << pop_history[40].N << std::endl;

    epoch_40_to_41 epoch_40(epoch_39, pop_history[41], inflection_points[41]);
    std::cout << "Finished 41 demographic events" << std::endl;

    epoch_41_to_42 epoch_41(epoch_40, pop_history[42], inflection_points[42]);
       std::cout << "Finished 42 demographic events" << std::endl;

    epoch_42_to_43 epoch_42(epoch_41, pop_history[43], inflection_points[43]);
        std::cout << "Finished 43 demographic events" << std::endl;

    epoch_43_to_44 epoch_43(epoch_42, pop_history[44], inflection_points[44]);
        std::cout << "Finished 44 demographic events" << std::endl;

    epoch_44_to_45 epoch_44(epoch_43, pop_history[45], inflection_points[45]);
        std::cout << "Finished 45 demographic events" << std::endl;
        std::cout << "Population size at epoch 45: " << pop_history[45].N << std::endl;

    epoch_45_to_46 epoch_45(epoch_44, pop_history[46], inflection_points[46]);
    std::cout << "Finished 46 demographic events" << std::endl;
    epoch_46_to_47 epoch_46(epoch_45, pop_history[47], inflection_points[47]);
    epoch_47_to_48 epoch_47(epoch_46, pop_history[48], inflection_points[48]);
    epoch_48_to_49 epoch_48(epoch_47, pop_history[49], inflection_points[49]);
    epoch_49_to_50 epoch_49(epoch_48, pop_history[50], inflection_points[50]);
    std::cout << "Finished 50 demographic events" << std::endl;
    epoch_50_to_51 epoch_50(epoch_49, pop_history[51], inflection_points[51]);
    epoch_51_to_52 epoch_51(epoch_50, pop_history[52], inflection_points[52]);
    epoch_52_to_53 epoch_52(epoch_51, pop_history[53], inflection_points[53]);
    epoch_53_to_54 epoch_53(epoch_52, pop_history[54], inflection_points[54]);
    epoch_54_to_55 epoch_54(epoch_53, pop_history[55], inflection_points[55]);
    epoch_55_to_56 epoch_55(epoch_54, pop_history[56], inflection_points[56]);

    // Migration parameters, no--migration
    mig_const mig_model;

    // Selection parameters
    Sim_Model::selection_constant weak_del((float)sel_coef);

    // SFS parameters
    int sample_size = num_samples; // number of samples in SFS
    int num_iter = 2;             // number of iterations
    Spectrum::SFS my_spectra;

    cudaEvent_t start, stop; // CUDA timing functions
    float elapsedTime;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    float avg_num_mutations = 0;
    float avg_num_mutations_sim = 0;
    std::vector<std::vector<float>> results(num_iter); // storage for SFS results
    for (int j = 0; j < num_iter; j++)
    {
        results[j].reserve(sample_size);
    }
    std::cout << "Starting iteration 1 of simulation" << std::endl;
    for (int j = 0; j < num_iter; j++)
    {
        if (j == num_iter / 2)
        {
            cudaEventRecord(start, 0);
        } // use 2nd half of the simulations to time simulation runs + SFS creation

        b.sim_input_constants.seed1 = 0xbeeff00d + 2 * j; // random number seeds
        b.sim_input_constants.seed2 = 0xdecafbad - 2 * j;
        GO_Fish::run_sim(b, mutation, epoch_55, mig_model, weak_del, outbred, codominant, Sim_Model::bool_off(), Sim_Model::bool_off());
        Spectrum::site_frequency_spectrum(my_spectra, b, 0, 0, sample_size);

        avg_num_mutations += ((float)my_spectra.num_mutations) / num_iter;
        avg_num_mutations_sim += b.maximal_num_mutations() / num_iter;

        for (int i = 0; i < sample_size; i++)
        {
            results[j][i] = my_spectra.frequency_spectrum[i];
        }
    }

    elapsedTime = 0;
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // output SFS simulation results
    std::cout << "SFS :" << std::endl
              << "allele count\tavg# mutations\tstandard dev\tcoeff of variation (aka relative standard deviation)" << std::endl;
    for (int i = 1; i < sample_size; i++)
    {
        double avg = 0;
        double std = 0;
        float num_mutations;
        for (int j = 0; j < num_iter; j++)
        {
            num_mutations = b.num_sites() - results[j][0];
            avg += results[j][i] / (num_iter * num_mutations);
        }
        for (int j = 0; j < num_iter; j++)
        {
            num_mutations = b.num_sites() - results[j][0];
            std += 1.0 / (num_iter - 1) * pow(results[j][i] / num_mutations - avg, 2);
        }
        std = sqrt(std);
        std::cout << i << "\t" << avg << "\t" << std << "\t" << (std / avg) << std::endl;
    }

    std::ofstream output_file("./out_sfs_57_epoch.txt");
    for (int i = 0; i < sample_size; i++)
    {
        output_file << my_spectra.frequency_spectrum[i] << "\n";
    }

    std::cout << "\nnumber of sites in simulation: " << b.num_sites() << "\ncompact interval: " << b.last_run_constants().compact_interval;
    std::cout << "\naverage number of mutations in simulation: " << avg_num_mutations_sim << "\naverage number of mutations in SFS: " << avg_num_mutations << "\ntime elapsed (ms): " << 2 * elapsedTime / num_iter << std::endl;
}

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

int main(int argc, char **argv)

{
    // this is the mutation rate scaled with respect to number of sites, mutation_rate*(number of sites)
    float mut_rate = 0.3426;
    // this is a point selection coefficient the selection coefficient will remain the same for the population, this is the un-scaled selection coefficient
    float PointSel = -.005;
    int num_samples = 113000;

    // Number of samples for to generate for the site-frequency spectrum (SFS

    // Eventually this will read in a demographic history file for easier command line use instead of having to re-compile for every new demography <- possible but will still require a compilation step as Functors (functions passed as templates) need to be known at compile time (a requirement of GPUs), I have not yet added the ability to do this to the library, technically there are other libraries that will allow this, but I haven't merged them with my API to make it easy. It's on my TODO list

    if (argc != 4) // 3 Total parameters, [executable, scaled mutation rate, unscaled selection coefficient, num_samples]
    {
        fprintf(stderr, "Error: The number of arguments given in the command line is not correct. In this version you need to pass in a selection cofficient and unscaled mutation rate, format is: ./GOFish scaled_mutation_rate unscaled_selection coefficient num_samples \n");
        // exit(8);
        std::cout << "Using default values" << std::endl;
    }
    else
    {

        mut_rate = atof(argv[1]);
        PointSel = atof(argv[2]);
        num_samples = atoi(argv[3]);
    }

    std::cout << "Scaled Mutation Rate: " << mut_rate << std::endl;
    std::cout << "Inscaled Point Selection: " << PointSel << std::endl;
    std::cout << "Number of samples to generate SFS: " << num_samples << std::endl;

    std::cout << "Running simulations" << std::endl;

    run_validation_test(mut_rate, PointSel, num_samples);
}
