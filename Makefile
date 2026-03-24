# LEGO - Lattice Elemental Geometry Operations
# Makefile for building the lego crystal structure generator

CC       ?= gcc
CFLAGS   ?= -O2 -Wall -Wextra
LDFLAGS  ?= -lm -lz

# OpenMP: auto-detect platform
# macOS clang doesn't support -fopenmp by default.
# Use: brew install libomp && make CC=gcc-14
# Or build without OpenMP: make OMP_FLAG=""
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    OMP_FLAG ?=
else
    OMP_FLAG ?= -fopenmp
endif

# Where built-in POSCAR files are installed
STRUCTURES_DIR ?= structures

SRCS     := $(wildcard src/*.c)
OBJS     := $(SRCS:.c=.o)
TARGET   := lego

# Pass structures directory as compile-time default
DEFINES  := -DSTRUCTURES_DIR=\"$(STRUCTURES_DIR)\"

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OMP_FLAG) -o $@ $^ $(LDFLAGS)

src/%.o: src/%.c src/lego.h
	$(CC) $(CFLAGS) $(OMP_FLAG) $(DEFINES) -c -o $@ $<

clean:
	rm -f src/*.o $(TARGET)

install: $(TARGET)
	@echo "Install not yet configured. Copy 'lego' binary and 'structures/' to your preferred location."
