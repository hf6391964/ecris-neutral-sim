LIBNAME = libneutrals
INCDIRS = -I../CGAL-4.9/include -I./include -I../gsl-2.2.1
CC = g++
CFLAGS_COMMON = -Wall -Wextra -std=c++11 -O2
CFLAGS = -g
#CFLAGS = 
PRECOMPILED_HEADER = include/cgal_and_typedefs.h
PRECOMPILED_HEADER_TARGET = $(PRECOMPILED_HEADER).gch

.PHONY: default all clean

default: $(LIBNAME)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard lib/*.cpp))
HEADERS = $(wildcard include/*.h)

$(PRECOMPILED_HEADER_TARGET): $(PRECOMPILED_HEADER)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS) $(INCDIRS) $(PRECOMPILED_HEADER) -o $(PRECOMPILED_HEADER_TARGET)

lib/%.o: lib/%.cpp $(HEADERS) $(PRECOMPILED_HEADER_TARGET)
	$(CC) $(CFLAGS_COMMON) $(CFLAGS) $(INCDIRS) -include $(PRECOMPILED_HEADER) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(LIBNAME): $(OBJECTS)
	ar cr $@.a $(OBJECTS)

clean:
	rm -f $(OBJECTS)
	rm -f $(LIBNAME).a
	rm -f $(PRECOMPILED_HEADER_TARGET)

