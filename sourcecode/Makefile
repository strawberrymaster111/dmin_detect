# My compiler:
CC = xlC

# Virtual path for include files:
#VPATH = ../include

# Libraries:
LIB = -lcomplex -L/zurich/usr/matlab4/curr/extern/lib/ibm_rs -lm

# Compiler flags:
# ---------------
# Normal setting:
#CCFLAGS = -I../include
#
# For debugging (Note that -g disables inlining):
CCFLAGS = -I../include -g -qinfo=nocnv:notrd -qcheck -bloadmap -bnoquiet

#
# For optimization:
#CCFLAGS = -I../include -O
.SUFFIXES: .C


OBJ = Main.o BlockCodeGenerator.o \
SoftDemodulator.o ChannelTransition.o ParityCheck.o FactorGraph.o BlaumSpectra.o Node.o FunctionNode.o VariableNode.o LocalFunction.o ParityCheckNode.o AWGN.o


# How to build the project:
# -------------------------
lib:	    Main 


#mytest:      mytest.o random.o codes.o spadet.o ldpc.o
#	$(CC) $(CCFLAGS) -o mytest mytest.o random.o codes.o spadet.o ldpc.o $(LIB)
#mytest.o:    mytest.C iofnct.h random.h codes.h spadet.h ldpc.h

Main: Main.o $(OBJ)
	$(CC) $(CFLAGS) -o Main $(OBJ) $(LIB) -bmaxdata:0x80000000



# Dependencies:
# -------------

Main.o: BlockCodeGenerator.h FactorGraph.h SoftDemodulator.h AWGN.h
	$(CC) $(CFLAGS) -c Main.C

BlockCodeGenerator.o: BlockCodeGenerator.h
	$(CC) $(CFLAGS) -c BlockCodeGenerator.C

BlaumSpectra.o: BlaumSpectra.h
	$(CC) $(CFLAGS) -c BlaumSpectra.C

SoftDemodulator.o: SoftDemodulator.h Constant.h
	$(CC) $(CFLAGS) -c SoftDemodulator.C

ChannelTransition.o: ChannelTransition.h FunctionNode.o
	$(CC) $(CFLAGS) -c ChannelTransition.C

ParityCheck.o: ParityCheck.h FunctionNode.o
	$(CC) $(CFLAGS) -c ParityCheck.C

FactorGraph.o: FactorGraph.h ParityCheckNode.o ChannelTransition.o VariableNode.o LocalFunction.o ParityCheck.o
	$(CC) $(CFLAGS) -c FactorGraph.C

Node.o: Node.h
	$(CC) $(CFLAGS) -c Node.C

FunctionNode.o: FunctionNode.h Constant.h Node.o
	$(CC) $(CFLAGS) -c FunctionNode.C

VariableNode.o: VariableNode.h Constant.h Node.o
	$(CC) $(CFLAGS) -c VariableNode.C

LocalFunction.o: LocalFunction.h 
	$(CC) $(CFLAGS) -c LocalFunction.C

ParityCheckNode.o: ParityCheckNode.h FunctionNode.o
	$(CC) $(CFLAGS) -c ParityCheckNode.C

AWGN.o: AWGN.h 
	$(CC) $(CFLAGS) -c AWGN.C


# Clean facility:
#---------------
clean:     clean_objects

clean_objects:     `rm *.o`











