CC = gcc
CC = mpicc
CFLAGS =-c -w -Wall -Wextra -pedantic -pg
LDFLAGS = -lgsl -lgslcblas -lm
EXECUTABLE = main_brownconf
head_inc = dummy

#$(head_inc): header_gen.py
#	python header_gen.py $(INT) $(CONF)


$(EXECUTABLE): main_brownconf.o
	$(CC) -w -Wall  -Wextra -pedantic main_brownconf.o par_sim.c conf_$(CONF).c int_$(INT).c $(LDFLAGS) -o $(EXECUTABLE)
	rm *.o

main_brownconf.o: main_brownconf.c
	$(CC) $(CFLAGS) main_brownconf.c


