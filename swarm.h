#pragma once

#include<iostream>
#include<cstdlib>
#include<vector>
#include<limits>
#include<cmath>

/*
	Each particle is a culmination of the weights used for the neural network
  */
class Swarm{
public:
	class particle{
	public:
		double fitness; // more fit if it is small
		std::vector<double> weights;
		std::vector<double> best_pweight;
		std::vector<double> velocity;
		int network[6] = {1, 5, 25, 25, 5, 1};
		//int network[4] = {1, 3, 3, 1};
		std::vector<std::vector<std::vector<double>>> layer_list;
		
		double run_network(double input);

		particle();
	};

	std::vector<particle> swarm_vec;
	particle best_gweight;
	double best_gweight_fitness;

	int iterations;
	int particle_amt;
	double w; // inertia weight
	double w_delta = w/iterations;
	double c1; // influence of personal best position
	double c2; // influence of global best position
	
	//these function update every particle
	void fitness_function(particle* p);
	void update_velocity(particle* p);
	void update_weight(particle* p);

	double siner(double input){
		return best_gweight.run_network(input);
	}

	void iterate_swarm(); // this is where we weight update

	Swarm(int i, int p, double weight, double c1_, double c2_);
};
