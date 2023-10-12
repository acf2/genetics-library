TESTDIR=test
BINDIR=bin

TESTSRC := $(shell ls $(TESTDIR)/*.cpp)
TESTEXEC := $(patsubst $(TESTDIR)/%.cpp,$(BINDIR)/%,$(TESTSRC))

CPPFLAGS=-std=c++20
DEBUG=-DDEBUG -g -w -pedantic -Wall
CPPFLAGS=-DNDEBUG -O2

INC=-I./include

all: CPPFLAGS := $(CPPFLAGS) $(RELEASE)
all: test

debug: CPPFLAGS := $(CPPFLAGS) $(DEBUG)
debug: test

clean:
	-rm -r $(BINDIR)

test: $(TESTEXEC)

$(BINDIR)/%: $(TESTDIR)/%.cpp
	g++ $(CPPFLAGS) $(INC) -o $@ $< 

$(TESTEXEC) : | $(BINDIR)

$(BINDIR):
	mkdir $(BINDIR)

.PHONY: all clean install test
