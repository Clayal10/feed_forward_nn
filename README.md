# feed_forward_nn
This is a feed forward neural network I created to approximate the output of the sine function. I am using Particle Swarm Optimization as well.
## Building
There is a makefile included to compile with g++. I am using an O3 optimizer to make sure it won't take too long ensuring enough iterations are ran.
## Other applications
It would be easy enough to change what you're using this for. You would need to change the inputs and what function is being estimated. I am using the built in cmath sin() function.
## Output
A CSV file is the output. It includes a column for the input, the network output, and the real sin() output.
