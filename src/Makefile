OBJS = main.o swarm.o
EXFLAGS = -O3 -o run

try: ${OBJS}
	g++ ${OBJS} ${EXFLAGS} -g

clean:
	rm -f run ${OBJS}

test: all
	./run


