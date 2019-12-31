#*************************************************************************#
#                                                                         #
#  Program: Makefile                                                      #
#  Version: 2.0                                                           #
#  By: Mette Olufsen                                                      #
#  Date: 14. Jan. 1997                                                    #
#  Last updated by: MJ Colebank, Dec. 31, 2019                            #
#                                                                         # 
#  A makefile that ensures that all modules are linked together in the    #
#  right order.                                                           #
#*************************************************************************#

CXX=g++
CXXFLAGS=-O2 -g -Wall -D_REENTRANT


LIBS=$(FLIBS) -lm

LDFLAGS=-O2

OBJS=tools.o sor06.o arteries.o

MAIN=sor06

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CXX) -o $(MAIN) $(LDFLAGS) $(OBJS) $(LIBS)
	
sor06.o: sor06.c sor06.h
	$(CXX) -c $(CXXFLAGS) sor06.c
	
arteries.o: arteries.c arteries.h tools.h sor06.h
	$(CXX) -c $(CXXFLAGS) arteries.c
	
tools.o: tools.c tools.h
	$(CXX) -c $(CXXFLAGS) tools.c		
	
clean:
	-rm -f *.o *.mod Zhat* *.2d 
	
veryclean: clean
	-rm $(MAIN) a.out *~ sor06
