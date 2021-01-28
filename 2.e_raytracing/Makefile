main: main_all.o Global.o precalculate.o Trace_rice.o
	gcc -g -Wall -lm main_all.o Global.o precalculate.o Trace_rice.o -o test

main_all.o: main_all.c Defs.h Global.h Trace.h precalculate.h
	gcc -g -Wall -c main_all.c

Global.o: Global.c Defs.h
	gcc -g -Wall -c Global.c

precalculate.o: precalculate.c Defs.h Global.h Trace.h
	gcc -g -Wall -c precalculate.c

Trace_rice.o: Trace_rice.c Defs.h Global.h
	gcc -g -Wall -c Trace_rice.c

clean:
	@echo "cleanning project"
	-rm *.o
	@echo "complete"
