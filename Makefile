TARGET = sim
INCDIRS = -I../CGAL-4.9/include #-I../gsl-2.2.1
LIBS = -lpthread
STATIC_LIBS = ../CGAL-4.9/lib/libCGAL.a #../gsl-2.2.1/.libs/libgsl.a
CC = g++
CFLAGS_COMMON = -Wall -std=c++11 -O3 #-flto
# CFLAGS = -g
CFLAGS = #-O3

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS) $(INCDIRS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS_COMMON) -o $@ $(LIBS) $(OBJECTS) $(STATIC_LIBS) 

clean:
	rm -f *.o
	rm -f $(TARGET)

