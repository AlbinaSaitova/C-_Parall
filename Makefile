CFLAGS=-O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Wcast-qual -Wno-suggest-attribute=format -pthread
all: a.out
a.out: task02.o solve.o HW2.h
	g++ $(CFLAGS) $^ -o $@
task02.o: task02.cpp  HW2.h
	g++ -c $(CFLAGS) $<
solve.o:  solve.cpp HW2.h
	g++ -c $(CFLAGS) $<
clean:
	rm -f *.out *.o *.bak *.txt	
