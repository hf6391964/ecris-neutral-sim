TARGET = sim
INCDIRS = -I../CGAL-4.8.1/include -I../gsl-2.2.1
LIBS = 
STATIC_LIBS = ../CGAL-4.8.1/lib/libCGAL.a ../gsl-2.2.1/.libs/libgsl.a
CC = g++
CFLAGS_COMMON = -Wall -std=c++11
# CFLAGS = -g
CFLAGS = -O2

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS) $(INCDIRS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(STATIC_LIBS) -Wall $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f $(TARGET)

