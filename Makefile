# Compile flags
# ====================================================================
# Warning flags
flag_w01 = -Wall -Wextra -Wconversion -pedantic -Wno-unused-parameter
# Debugging flags
flag_d01 = -g
# Optimization flags
flag_o01 = -O3 -ftree-vectorize -ftree-loop-vectorize -funroll-loops
flag_o02 = -march=native -Ofast -ffast-math
flag_o03 = -fopt-info-vec -ftree-vectorizer-verbose=2

CFLAGS = ${flag_o01} ${flag_o02} ${flag_w01} ${flag_d01} ${flag_o03}

# Directories
SRCDIR = src
INCDIR = include
BUILDDIR = build

# Main program (to be set via command line or manually)
MAIN = 

# Libraries
LDLIBS = -lgsl -lgslcblas -lm

# Header files (space-separated list)
HEADERS = parameters ini common

# Transform header list into .o and .h file lists
comma := ,
HEADERS_LIST := $(subst  ,$(comma),$(HEADERS))
HEADERS_LIST := $(subst $(comma), ,$(HEADERS_LIST))
HEADERS_O := $(addsuffix .o,$(HEADERS_LIST))
HEADERS_H := $(addsuffix .h,$(HEADERS_LIST))

# Object files (for dependencies)
OBJS = $(HEADERS_O:%=$(BUILDDIR)/%)

# Rules
.PHONY: all clean

all: $(BUILDDIR)/$(MAIN)

# Rule to build the main program
$(BUILDDIR)/$(MAIN): $(SRCDIR)/$(MAIN).o $(OBJS)
	@mkdir -p $(BUILDDIR)
	@$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

# Rule to compile source files into object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HEADERS_H:%=$(INCDIR)/%)
	@mkdir -p $(BUILDDIR)
	@$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

# Clean up
clean:
	@rm -f $(BUILDDIR)/* $(SRCDIR)/*.o

.PHONY: clean all
