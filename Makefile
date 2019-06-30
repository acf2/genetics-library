LIBNAME=gntc
TESTDIR=test
BINDIR=bin

TESTSRC := $(shell ls $(TESTDIR)/*.cpp)
TESTEXEC := $(patsubst $(TESTDIR)/%.cpp,$(BINDIR)/%,$(TESTSRC))

#CPPFLAGS=-std=c++14 -DDEBUG -g -w -pedantic -Wall
CPPFLAGS=-std=c++14 -DNDEBUG -O2

INC=-I./include

all: test

clean:
	-rm -r $(BINDIR)

install:
	cp include/* /home/acf2/0000/gear/NSL/include/

test: $(TESTEXEC)

$(BINDIR)/%: $(TESTDIR)/%.cpp
	g++ $(CPPFLAGS) $(INC) -o $@ $< 

$(TESTEXEC) : | $(BINDIR)

$(BINDIR):
	mkdir $(BINDIR)

.PHONY: all clean install test
