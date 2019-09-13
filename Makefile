#*************************************************************************#
#                                                                         #
#  Program: Makefile                                                      #
#  Version: 2.0                                                           #
#  By: Mette Olufsen                                                      #
#  Date: 14. Jan. 1997                                                    #
#                                                                         # 
#  A makefile that ensures that all modules are linked together in the    #
#  right order.                                                           #
#*************************************************************************#

# $Id: Makefile,v 1.8 2010-10-20 14:39:59 heine Exp $

CXX=g++
CXXFLAGS=-O2 -g -Wall -D_REENTRANT


LIBS=$(FLIBS) -lm

LDFLAGS=-O2

OBJS=tools.o sor06.o arteries.o

MAIN=sor06

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CXX) -o $(MAIN) $(LDFLAGS) $(OBJS) $(LIBS)
	
sor06.o: sor06.C sor06.h
	$(CXX) -c $(CXXFLAGS) sor06.C
	
arteries.o: arteries.C arteries.h tools.h sor06.h
	$(CXX) -c $(CXXFLAGS) arteries.C
	
tools.o: tools.C tools.h
	$(CXX) -c $(CXXFLAGS) tools.C		
	
clean:
	-rm -f *.o *.mod Zhat* *.2d
	
veryclean: clean
	-rm $(MAIN) a.out *~ sor06
