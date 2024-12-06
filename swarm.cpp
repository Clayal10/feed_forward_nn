#include "swarm.h"

double rand_double(double min, double max){
	return min+(max - min)*(rand()%RAND_MAX) / RAND_MAX;
}

/*activation functions*/
double unipolar(double x){
	return 1 / (1+pow(M_E, -x));
}
double bipolar(double x){
	return (2 / (1+pow(M_E, -x))) - 1;
}
double relu(double x){
	if(x >=0) return x;
	return 0;
}

/*helper*/
void print_matrix(std::vector<std::vector<double>> mat){
	for(int i=0; i<mat.size(); i++){
		for(int j=0; j<mat[0].size(); j++){
			std::cout << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void transpose(std::vector<std::vector<double>>* layer){
	std::vector<double> buffer;
	std::vector<std::vector<double>> mat_buffer;
	//transpose layer
	for(int i=0; i<layer->at(0).size(); i++){
		for(int j=0; j<layer->size(); j++){
			buffer.push_back(layer->at(j)[i]);
		}
		mat_buffer.push_back(buffer);
		buffer.clear();
	}

	
	layer->clear();

	for(int i=0; i<mat_buffer.size(); i++){
		for(int j=0; j<mat_buffer[0].size(); j++){
			buffer.push_back(mat_buffer[i][j]);
		}
		layer->push_back(buffer);
		buffer.clear();
	}

}


/*Matrix Multiplication*/
//will need to be careful about orientation;
void matrix_mul(std::vector<std::vector<double>>* in, std::vector<std::vector<double>> layer){
	/*translate 2nd matrix*/
	std::vector<double> buffer;
	std::vector<std::vector<double>> mat_buffer;

	//conditional for transpose
	transpose(in);
	if(in->at(0).size() != layer.size()){
		transpose(&layer);
	}
	//std::cout << in->size() << "x" << in->at(0).size() << "\t*\t" << layer.size() << "x" << layer[0].size() << std::endl; 
	
	buffer.clear(); // need to use this again
	mat_buffer.clear();
	double buff = 0;
	for(int k=0; k<layer[0].size(); k++){
		for(int i=0; i<in->size(); i++){
			for(int j=0; j<in->at(0).size(); j++){
				buff += in->at(i)[j] * layer[j][k];
			}
			buffer.push_back(buff);
			buff = 0;
		}
		mat_buffer.push_back(buffer); // be 1 vector with 5 element of 5 vector with 1 element each
		buffer.clear();
	}
	in->clear();

	for(int i=0; i<mat_buffer.size(); i++){
		for(int j=0; j<mat_buffer[0].size(); j++){
			buffer.push_back(mat_buffer[i][j]);
		}
		in->push_back(buffer);
		buffer.clear();
	}
	// 'in' should be updated with the old input * weight matrix layer
	

}

void apply_bias(std::vector<double>* v, std::vector<double> bias){
	for(int i=0; i<v->size(); i++){
		v->at(i) += bias[i];
	}
}

/*Class methods*/
Swarm::Swarm(int i, int p, double weight, double c1_, double c2_){
	iterations = i;
	particle_amt = p;
	w = weight;
	c1 = c1_;
	c2 = c2_;
	best_gweight_fitness = std::numeric_limits<double>::max();

	for(int i = 0; i<particle_amt; i++){
		swarm_vec.push_back(particle());
	}
	
	for(int i=0; i<swarm_vec.size(); i++){
		fitness_function(&swarm_vec[i]);
	}

	
}

/*Updating velocity and weights*/
void Swarm::update_velocity(particle* p){ // will be ran every iteration, updates all particles
	double new_velocity;
	double r1 = rand_double(0, 1);
	double r2 = rand_double(0, 1);

	for(int i=0; i<p->velocity.size(); i++){
		new_velocity = w * p->velocity[i] + c1 * r1 * (p->best_pweight[i] - p->weights[i]) + c2 * r2 * (best_gweight.weights[i] - p->weights[i]);
		
		//bounds check
		if(new_velocity < -1) new_velocity = -1;
		else if(new_velocity > 1) new_velocity = 1;
		//set
		//std::cout << "old: " << p->velocity[i] << "\tnew: " << new_velocity << std::endl;

		p->velocity[i] = new_velocity;
	}

}

void Swarm::update_weight(particle* p){
	double new_weight;
	int cnt = 0;
	for(int i=0; i<p->layer_list.size(); i++){
		for(int j=0; j<p->layer_list[i].size(); j++){
			for(int k=0; k<p->layer_list[i][j].size(); k++){
				new_weight = p->layer_list[i][j][k] + p->velocity[cnt];
				//could check bounds?
				p->weights[cnt] = new_weight;
				p->layer_list[i][j][k] = new_weight;
				cnt++;
			}
		}
	}

}


Swarm::particle::particle(){
	fitness = std::numeric_limits<double>::max();

	int weight_amt = 0;
	for(int i=0; i<6; i++){
	//for(int i=0; i<5; i++){
	//for(int i=0; i<4; i++){
		weight_amt += network[i];
	}
	for(int i=0; i<weight_amt; i++){
		weights.push_back(rand_double(-1.0, 1.0));
		velocity.push_back(rand_double(-1.0, 1.0));
	}
	std::vector<std::vector<double>> matrix_buffer;
	std::vector<double> buffer;
	int j=0;
	for(int i=0; i<6; i++){
	//for(int i=0; i<5; i++){
	//for(int i=0; i<4; i++){
		
		if(i != 2 && i!= 3){
			for(int k=0; k<network[i]; k++){
				buffer.push_back(weights[j]);
				//std::cout << weights[j] << std::endl;
				j++;
			}
			matrix_buffer.push_back(buffer);
			buffer.clear();
		}else{
			for(int k=0; k<5; k++){
				for(int l=0; l<5; l++){
					buffer.push_back(weights[j]);
					//std::cout << weights[j] << std::endl;
					j++;
				}
				matrix_buffer.push_back(buffer);
				buffer.clear();
			}
		}
		layer_list.push_back(matrix_buffer);
		matrix_buffer.clear();

	}
	
}

void Swarm::fitness_function(particle* p){
	//calculate MSE
	double error_buff=0;
	double real, predicted;
	int counter = 0;
	for(double i = -3*M_PI; i<3*M_PI; i+= 0.05){
		predicted = p->run_network(i);
		
		real = sin(i);
		error_buff += (predicted - real) * (predicted - real);
		counter++;	
	}
	error_buff = error_buff/counter;
	if(error_buff < p->fitness){
		p->best_pweight = p->weights;
	}
	p->fitness = error_buff;

	if(p->fitness < best_gweight_fitness){
		best_gweight_fitness = p->fitness;
		best_gweight = *p;
	}

}

double Swarm::particle::run_network(double input){
	//this is where the input is ran through the network
	std::vector<std::vector<double>> in;
	in.push_back({input}); // just a single element for input

	for(int i=0; i<layer_list.size(); i++){ // this will loop through each layer, output being new matrix
		matrix_mul(&in, layer_list[i]);
		if(i != layer_list.size()-1){
			for(int j=0; j<in.size(); j++){
				for(int k=0; k<in[0].size(); k++){
					in[j][k] = bipolar(in[j][k]); // run each element through activation function
				}
			}
		}
	}
	return in[0][0];
}

/*The Iteration*/
void Swarm::iterate_swarm(){ // at this point, the swarm is initialized
	int count = 0;
	while(count < iterations){
		for(int i=0; i<swarm_vec.size(); i++){// run through every particle
			update_velocity(&swarm_vec[i]);
			update_weight(&swarm_vec[i]); // weight is updated here
			fitness_function(&swarm_vec[i]);

		}
		//if(w - w_delta > 0)
		//	w -= w_delta/2;
		//this is used for printing the weights before and after training
		if(count == 0 || count == iterations-1){
			for(int i=0; i<best_gweight.layer_list.size(); i++){
				print_matrix(swarm_vec[0].layer_list[i]);
				std::cout << "\n";
			}
			std::cout << "\n";
		}
		count++;


		if(count%10 == 0){
			std::cout << count << std::endl;
		}
	}
}

