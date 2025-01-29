# feed_forward_nn
This is a feed forward neural network I created to approximate the output of the sine function. I am using Particle Swarm Optimization as well.
## Building
Build with the included Makefile. Uses O3 optimizer to speed it up a bit.
## Other applications
It would be easy enough to change what you're using this for. You would need to change the inputs and what function is being estimated. I am using the built in cmath sin() function.
## Output
A CSV file is the output. It includes a column for the input, the network output, and the real sin() output.
