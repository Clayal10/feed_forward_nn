try:
	g++ main.cpp swarm.cpp -O3 -g
test: all
	a.out
