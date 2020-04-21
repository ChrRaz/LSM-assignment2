# Makefile
#

TARGET = poisson

OBJECTS	= print.o init.o

# options and settings for the GCC compilers
#
CC			 = mpicc
DEFS		 =
OPT			 = -g -O3
IPO			 =
ISA			 =
CHIP		 =
ARCH		 =
PARA		 = -fopenmp
CFLAGS	 = $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS)
LDFLAGS	 = -lm
# OMPFLAGS = -fopenmp

TEST_SIZE       = 100
TEST_ITER       = 15000
TEST_THRESHOLD  = 0.01
TEST_START_T    = 0
TEST_OUTPUT     = 4
TESTFLAGS       = $(TEST_SIZE) $(TEST_ITER) $(TEST_THRESHOLD) $(TEST_START_T) $(TEST_OUTPUT)

.PHONY: all clean realclean bench cleanbench
all: $(TARGET)

poisson: main.o $(OBJECTS) jacobi.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

main.o:
	$(CC) -o $@ $(CFLAGS) -c main.c

clean:
	@/bin/rm -vf core *.o *~

realclean: clean
	@/bin/rm -vf $(TARGET_J)  $(TARGET_GS)

# DO NOT DELETE

main.o: main.c print.h jacobi.h
print.o: print.h
init.o: init.h
