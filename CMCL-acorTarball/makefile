acc.o: acc.cpp acor.h makefile
	g++ -c acc.cpp

acrc.o: acrc.cpp acor.h makefile
	g++ -c acrc.cpp

acor: acc.o acrc.o acor.h
	g++ -o acor acc.o acrc.o 
   
acorTest: acorTest.cpp makefile            # the tester is its own thing.
	g++ -o acorTest acorTest.cpp 

acorTarball: makefile README AcorWriteup.tex AcorWriteup.pdf acc.cpp acor.h acrc.cpp acorTest.cpp 
	tar -cf acorTarball.tar makefile README AcorWriteup.tex AcorWriteup.pdf acc.cpp acor.h acrc.cpp acorTest.cpp