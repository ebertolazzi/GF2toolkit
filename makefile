SRCS = \
srcs/GF2toolkit_4RussianMult10bit.cc \
srcs/GF2toolkit_4RussianMult11bit.cc \
srcs/GF2toolkit_4RussianMult2bit.cc \
srcs/GF2toolkit_4RussianMult3bit.cc \
srcs/GF2toolkit_4RussianMult4bit.cc \
srcs/GF2toolkit_4RussianMult5bit.cc \
srcs/GF2toolkit_4RussianMult6bit.cc \
srcs/GF2toolkit_4RussianMult7bit.cc \
srcs/GF2toolkit_4RussianMult8bit.cc \
srcs/GF2toolkit_4RussianMult9bit.cc \
srcs/GF2toolkit_4RussianMultClasses.cc \
srcs/GF2toolkit_4RussianMultReduced.cc \
srcs/GF2toolkit_4RussianMultReducedShift.cc \
srcs/GF2toolkit_4RussianTables.cc \
srcs/GF2toolkit_CompactColumn.cc \
srcs/GF2toolkit_InvertBlock.cc \
srcs/GF2toolkit_LU.cc \
srcs/GF2toolkit_LUcompute.cc \
srcs/GF2toolkit_LUcomputeRecurr.cc \
srcs/GF2toolkit_LUextract.cc \
srcs/GF2toolkit_Matrix.cc \
srcs/GF2toolkit_MatrixLinvApply.cc \
srcs/GF2toolkit_MatrixLower.cc \
srcs/GF2toolkit_MatrixMult.cc \
srcs/GF2toolkit_MatrixUinvApply.cc \
srcs/GF2toolkit_Strassen.cc \
srcs/GF2toolkit_m4ri.cc \
srcs/GF2toolkit_popCount.cc \
srcs/Random.cc \
srcs/TimeMeter.cc

OBJS = $(SRCS:.cc=.o)

DEPS = \
srcs/GF2toolkit_4Russian.hh \
srcs/GF2toolkit_4RussianMacros.hh \
srcs/GF2toolkit_InvertBlock.hh \
srcs/GF2toolkit_LU.hh \
srcs/GF2toolkit_Matrix.hh \
srcs/GF2toolkit_MatrixLower.hh \
srcs/GF2toolkit_MatrixMult.hh \
srcs/GF2toolkit_Strassen.hh \
srcs/GF2toolkit_common.hh \
srcs/GF2toolkit_m4ri.hh \
srcs/GF2toolkit_popCounts.hh \
srcs/Random.hh \
srcs/TimeMeter.hh

#
# SELECT COMPILERS AND FLAGS
#
#CC     = llvm-gcc
#CXX    = llvm-g++
#CC     = clang
#CXX    = clang++
CC     = gcc
CXX    = g++

INC    = -Isrcs -I/usr/local/include
CFLAGS = -Wall -O2 -funroll-loops -msse3 -msse2 -msse -mmmx -m64 -fPIC 
LIBS   = -L. -lGF2toolkit -lm4ri

#AR     = ar rcs
AR     = libtool -static -o 

all: libGF2toolkit.a
	$(CXX) $(CFLAGS) $(INC) -o test_n1  tests/test_n1_BitsOP.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n2  tests/test_n2_PopCount.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n3  tests/test_n3_baseOP.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n4  tests/test_n4_MatrixMatrixMult.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n4b tests/test_n4_testLib.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n5  tests/test_n5_timeMatrixMult.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n6  tests/test_n6_MatrixLower.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n7  tests/test_n7_InvertBlock.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n8  tests/test_n8_timeLinv.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n9  tests/test_n9_LUdecomposition.cc $(LIBS)
	$(CXX) $(CFLAGS) $(INC) -o test_n10 tests/test_n10_RankComputationTimeComparison.cc $(LIBS)

srcs/%.o: srcs/%.cc $(DEPS)
	$(CXX) $(CFLAGS) $(INC) -c $< -o $@ 

srcs/%.o: srcs/%.c $(DEPS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

libGF2toolkit.a: $(OBJS)
	$(AR) libGF2toolkit.a $(OBJS) 

doc:
	doxygen
	
clean:
	rm -f libGF2toolkit.a srcs/*.o test_n*