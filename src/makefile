CC = gcc
CFLAGS = -c -w -Wall -Wextra -pedantic
LDFLAGS = -lgsl -lgslcblas -lm
EXECUTABLE = main_brownconf

CONF = conf_splitter.c
INT = int_lennardjones.c

head_inc = dummy

#$(head_inc): header_gen.py
#	python header_gen.py $(INT) $(CONF)


$(EXECUTABLE): main_brownconf.o
	$(CC) -w -Wall  -Wextra -pedantic main_brownconf.o sim_config.c results_transport.c code_handling.c random_numb_gen.c simulation_core.c print_routines.c $(CONF) $(INT) $(LDFLAGS) -o $(EXECUTABLE)
	rm *.o

main_brownconf.o: main_brownconf.c
	$(CC) $(CFLAGS) main_brownconf.c
