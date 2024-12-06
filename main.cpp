#include<iostream>
#include<vector>
#include<cmath>
#include<cstdlib>
#include<fstream>

#include "swarm.h"

using namespace std;



int main(){
	srand((unsigned int)time(0));
	
	// Iter, size, inertia, c1, c2
	Swarm swarm(2000, 30, 0.7, 2, 2);
	
	cout << "Results and weights of the best particle before training the network:\n" <<
		"======================================================================\n";
	cout << "PI/2: " << swarm.siner(M_PI/2) << std::endl;
	cout << "PI: " << swarm.siner(M_PI) << std::endl;
	cout << "0: " << swarm.siner(0) << std::endl;
	cout << "best fitness: " << swarm.best_gweight.fitness << endl;
	swarm.iterate_swarm();
	cout << "Results and weights of the best particle after training the network:\n" <<
		"======================================================================\n";
	cout << "3PI/2: " << swarm.siner(3*M_PI/2) << std::endl;
	cout << "PI/2: " << swarm.siner(M_PI/2) << std::endl;
	cout << "PI: " << swarm.siner(M_PI) << std::endl;
	cout << "0: " << swarm.siner(0) << std::endl;
	cout << "best fitness: " << swarm.best_gweight.fitness << endl;


	ofstream fp;
	fp.open("data_NN.csv");
	
	fp << "Domain,Predicted,Actual\n";
	for(double i=M_PI*-3; i <= M_PI*3; i += 0.05){
		fp << i << "," << swarm.siner(i) << "," << sin(i) << "\n";
	}


	return 0;	
}


