CC = gcc
CFLAGS =-c
LDFLAGS = -lgsl -lgslcblas -lm
EXECUTABLE = main_brownconf


$(EXECUTABLE): main_brownconf.o
	$(CC) -w -Wall  -Wextra main_brownconf.o par_sim.c conf_$(CONF).c int_$(INT).c $(LDFLAGS) -o $(EXECUTABLE)
	rm *.o

main_brownconf.o: main_brownconf.c
	$(CC) $(CFLAGS) main_brownconf.c

