all: main inner_check trace_recal

main: main_all.o Global.o precalculate.o Trace_rice2.o
	gcc -g -Wall main_all.o Global.o precalculate.o Trace_rice2.o -o test -lm

main_all.o: main_all.c Defs.h Global.h Trace.h precalculate.h
	gcc -g -Wall -c main_all.c

inner_check: inner_check.o Global.o precalculate.o Trace_rice2.o
	gcc -g -Wall inner_check.o Global.o precalculate.o Trace_rice2.o -o inner_check -lm

inner_check.o: inner_check.c Defs.h Global.h Trace.h precalculate.h
	gcc -g -Wall -c inner_check.c

trace_recal: trace_recal.o Global.o precalculate.o
	gcc -g -Wall trace_recal.o Global.o precalculate.o -o trace_recal -lm

trace_recal.o: trace_recal.c Defs.h Global.h precalculate.h
	gcc -g -Wall -c trace_recal.c

Global.o: Global.c Defs.h
	gcc -g -Wall -c Global.c

precalculate.o: precalculate.c Defs.h Global.h Trace.h
	gcc -g -Wall -c precalculate.c

Trace_rice2.o: Trace_rice2.c Defs.h Global.h
	gcc -g -Wall -c Trace_rice2.c

clean:
	@echo "cleanning project"
	-rm *.o
	@echo "complete"
