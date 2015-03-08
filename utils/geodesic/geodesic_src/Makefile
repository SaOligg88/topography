#=================================================================
# 
#=================================================================
DEBUG=1

CC=g++
CFLAGS= -Wall -fPIC
LDFLAGS= -L. -lgeodesic  

ifeq ($(DEBUG),1)
	CFLAGS += -g
else
	CFLAGS += -O4
	LDFLAGS += -s 
endif

DEPS = geodesic_algorithm_base.h geodesic_algorithm_exact_elements.h \
       geodesic_algorithm_subdivision.h geodesic_memory.h \
       geodesic_algorithm_dijkstra_alternative.h \
       geodesic_algorithm_exact.h geodesic_constants_and_simple_functions.h \
       geodesic_mesh_elements.h geodesic_algorithm_dijkstra.h \
       geodesic_algorithm_graph_base.h geodesic_matlab_api.h                      geodesic_mesh.h


%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: all clean


all: lib example0 example1

lib: libgeodesic.so

libgeodesic.so: geodesic_matlab_api.o
	$(CC) -shared -Wl,-soname,$@ -o $@ $< -lc

example0: example0.o 
	$(CC) $< -o $@ $(LDFLAGS)

example1: example1.o 
	$(CC) $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o *.so example1 example2
