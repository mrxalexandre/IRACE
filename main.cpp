#include <iostream>
#include <sstream>
#include <vector>
#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include "bossa_timer.h"
#include "ArgPack.hpp"
#include <fstream>

std::vector<std::vector<double> > points;

int main(int argc, char* argv[]) {


	ArgPack single_ap(argc, argv);

//ArgPack::ap().time
	unsigned n = 0;			// size of chromosomes
	const unsigned p = ArgPack::ap().population;		// size of population
	const double pe = ArgPack::ap().populationElite;		// fraction of population to be the elite-set
	const double pm = ArgPack::ap().populationMutants;		// fraction of population to be replaced by mutants
	const double rhoe = ArgPack::ap().rhoe;	// probability that offspring inherit an allele from elite parent
	const unsigned K = ArgPack::ap().K;		// number of independent populations
	const unsigned MAXT = ArgPack::ap().threads;	// number of threads for parallel decoding

	const double cutoff_time = ArgPack::ap().time;

	// Reading instance
	std::string s;
	std::string file = ArgPack::ap().inputFile;

	ifstream f(file);

	int n_points, dim;
	//dstd::getline(std::cin, s);
	//std::stringstream st(s);
	f >> n_points >> dim;

	points = std::vector<std::vector<double> > (n_points);
	for(int i=0; i<n_points; ++i)
		points[i] = std::vector<double> (dim);


	//for(int d=0;d<dim;++d) {
		//CLUSTER_L[d] = 1000;
		//CLUSTER_U[d] = -1000;
	//}

	for(int i=0;i<n_points;++i) {

		//std::getline(std::cin, s);
		//std::stringstream st(s);
		for(int d=0; d<dim; ++d) {
			f >> points[i][d];

		}
	}



	SampleDecoder decoder;			// initialize the decoder



	const long unsigned rngSeed = ArgPack::ap().rngSeed;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	n = n_points;

	// initialize the BRKGA-based heuristic


	BossaTimer timer;
	timer.start();

	double bestValue = -1;
	double timerToBest;
	bool verbose = ArgPack::ap().verbose;
	unsigned clusters = ArgPack::ap().clusters;



	decoder.set_k(clusters);
	BRKGA<SampleDecoder, MTRand> algo(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);



	if(verbose)
		cout << "Initializing populations\n";




	unsigned generation = 0;		// current generation
	const unsigned X_INTVL =  ArgPack::ap().exchangeBest;	// exchange best individuals at every 100 generations


	const unsigned X_NUMBER = ArgPack::ap().exchangeTop;	// exchange top 2 best
	const unsigned MAX_GENS = ArgPack::ap().generations;	// run for 1000 gens

	do {

		++generation;


		if(verbose)
			cout << "Envolving generation " << generation << " ";





		algo.evolve();	// evolve the population for one generation

		double best_temp = (-1)* algo.getBestFitness();


		if(best_temp > bestValue) {
			timerToBest = timer.getTime();
			bestValue = best_temp;

			cout << "New best " << bestValue << endl;
		}

		if((++generation) % X_INTVL == 0) { // VEJA se usa isso
			algo.exchangeElite(X_NUMBER);	// exchange top individuals
		}



	} while (generation < MAX_GENS and timer.getTime() < cutoff_time);
	timer.pause();

//até aqui valores para gráfico ^^

		cout << "Best solution " << bestValue << endl;
		std::cout << "Total time = " << timer.getTime() << std::endl;
		std::cout << "Time to Best ttb = " << timerToBest << std::endl;



	return 0;
}
