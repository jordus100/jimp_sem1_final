aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o

wielomiany: main.o splines.o points.o aproks_wielomianem.o gaus/libge.a
	$(CC) -o wielomiany main.o splines.o points.o aproks_wielomianem.o -L gaus -l ge -lm

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

aproks_wielomianem.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c -lm aproks_wielomianem.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta wielomiany
