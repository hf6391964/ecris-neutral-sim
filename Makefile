LIBNAME = libneutrals
INCDIRS = -I../CGAL-4.9/include -I./include #-I../gsl-2.2.1
CC = g++
CFLAGS_COMMON = -Wall -Wextra -std=c++11 -O3 #-flto
# CFLAGS = -g
CFLAGS = 

.PHONY: default all clean

default: $(LIBNAME)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard lib/*.cpp))
HEADERS = $(wildcard include/*.h)

lib/%.o: lib/%.cpp $(HEADERS)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS) $(INCDIRS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(LIBNAME): $(OBJECTS)
	ar cr $@.a $(OBJECTS)

clean:
	rm -f $(OBJECTS)
	rm -f $(LIBNAME).a

