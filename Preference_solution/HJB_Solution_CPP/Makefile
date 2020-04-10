# define variables
VPATH   = ./
HDRDIR  = ./include

LIBDIR  = 

CXX     =  g++-9 -std=c++17 -fopenmp -ffast-math -fopt-info-vec-optimized 
LD      =  g++-9 -std=c++17 -fopenmp -ffast-math -fopt-info-vec-optimized



# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I$(HDRDIR) -I$(HDRDIR)/eigen3 -O3

# link flags to be used 
LDFLAGS = -I$(HDRDIR) -I$(HDRDIR)/eigen3 -L. -L./$(LIBDIR) -O3

# libraries to be linked in
LIBS =  

# types of files we are going to construct rules for
.SUFFIXES: .cpp 

# rule for .c files
.cpp.o:
	$(CXX) $(CFLAGS) -o $*.o -c $*.cpp

# list of objects to be compiled
OBJS    = \
        src/HJB.o\
        src/CG.o\
        src/HJB_Driver.o

HJBSolver:$(OBJS)
	$(LD)  $(LDFLAGS) -o HJBSolver $(OBJS) $(LIBS)

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)


