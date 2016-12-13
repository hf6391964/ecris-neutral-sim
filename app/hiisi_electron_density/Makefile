TARGET = main
INCDIRS = -I../../../CGAL-4.9/include -I../../include -I../../../gsl-2.2.1
LIBS = -lpthread
STATIC_LIBS = ../../../CGAL-4.9/lib/libCGAL.a ../../libneutrals.a ../../../gsl-2.2.1/.libs/libgsl.a
CC = g++
CFLAGS_COMMON = -Wall -Wextra -std=c++11 -O3
# CFLAGS = -g
CFLAGS = 

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS) $(INCDIRS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(CFLAGS_COMMON) $(CFLAGS) $(OBJECTS) $(STATIC_LIBS) $(LIBS)

clean:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

