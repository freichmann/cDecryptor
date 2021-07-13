CC=g++

COPTS=-std=c++1y -Wall -Wno-unused

LOPTS=-pthread

OOPTS=-O3 -ffast-math -fexpensive-optimizations

FILES = ./src/cDecryptor.cpp \
	./src/Options.h \
	./src/GaussianNorm.h ./src/GaussianNorm.cpp \
	./src/Lock.h ./src/Lock.cpp \
	./src/NGram.h ./src/NGram.cpp \
	./src/RatedScore.h ./src/RatedScore.cpp \
	./src/Score.h ./src/Score.cpp

cDecryptor: Makefile $(FILES); $(CC) $(COPTS) $(OOPTS) $(LOPTS) -o cDecryptor $(FILES);

all: cDecryptor