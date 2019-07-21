CXXFLAGS ?= -O2 -g -Wall
LDFLAGS ?=

SNDFILE_CFLAGS = $(shell pkg-config sndfile --cflags)
SNDFILE_LIBS = $(shell pkg-config sndfile --libs)

BBD_SOURCES := ../bbd_filter.cc ../bbd_filter.h ../bbd_line.cc ../bbd_line.h

all: bbd

clean:
	rm -f bbd

bbd: bbd.cc $(BBD_SOURCES)
	$(CXX) -std=c++11 $(CXXFLAGS) $(SNDFILE_CFLAGS) -o $@ $< $(LDFLAGS) $(SNDFILE_LIBS)

.PHONY: all clean
