#---Begin configuration---#
#---Begin configuration---#

prefix=/usr/local
GT_INSTALL_DIR=$(prefix)
GT_COMPILE_DIR=$(GT_INSTALL_DIR)/src/genometools

#----End configuration----#
#----End configuration----#

# Binaries
PE_EXE=bin/parseval
CN_EXE=bin/canon-gff3
VN_EXE=bin/vang
BINS=$(PE_EXE) $(CN_EXE) $(VN_EXE)

#----- Source, header, and object files -----#

# AEGeAn core class and module files
AGN_SRCS=$(wildcard src/core/Agn*.c)
AGN_OBJS=$(patsubst src/core/%.c,obj/%.o,$(AGN_SRCS))
AGN_HDRS=$(patsubst src/core%.c,inc/core/%.h,$(AGN_SRCS))

# ParsEval class and module files
PE_SRCS=$(wildcard src/ParsEval/Pe*.c)
PE_OBJS=$(patsubst src/ParsEval/%.c,obj/%.o,$(PE_SRCS))
PE_HDRS=$(patsubst src/ParsEval/%.c,inc/ParsEval/%.h,$(PE_SRCS))

# VAnG class and module files
VN_SRCS=$(wildcard src/VAnG/Vang*.c)
VN_OBJS=$(patsubst src/VAnG/%.c,obj/%.o,$(VN_SRCS))
VN_HDRS=$(patsubst src/VAnG/%.c,inc/VAnG/%.h,$(VN_SRCS))

# All class and module files
CLSS_MDL_OBJS=$(AGN_OBJS) $(PE_OBJS) $(VN_OBJS)

# Compilation settings
CC=gcc
CFLAGS=-Wall -DAGN_DATA_PATH='"$(prefix)/share/aegean"' # -Wno-unused-result
ifeq ($(cairo),no)
  CFLAGS += -DWITHOUT_CAIRO
endif
ifneq ($(errorcheck),no)
  CFLAGS += -Werror
endif
ifeq ($(optimize),yes)
  CFLAGS += -O3
endif
ifeq ($(64bit),yes)
  CFLAGS += -m64
endif
LDFLAGS=-lgenometools -lm -L$(GT_INSTALL_DIR)/lib
ifdef lib
  LDFLAGS += -L$(lib)
endif
INCS=-I $(GT_INSTALL_DIR)/include/genometools/ -I $(GT_COMPILE_DIR)/src -I inc/core -I inc/ParsEval -I inc/VAnG -I /usr/include/cairo/ -I /sw/include/cairo/

# Targets
all:		$(BINS)
		

install:	all
		@- test -d $(prefix)/bin || mkdir $(prefix)/bin
		cp $(BINS) $(prefix)/bin/.
		@- test -d $(prefix)/share || mkdir $(prefix)/share
		@- test -d $(prefix)/share/aegean || mkdir $(prefix)/share/aegean
		cp data/share/* $(prefix)/share/aegean/.

uninstall:	
		rm -r $(prefix)/$(PE_EXE)
		rm -r $(prefix)/share/parseval

clean:		
		rm -f $(BINS) $(CLSS_MDL_OBJS)

$(AGN_OBJS):	obj/%.o : src/core/%.c inc/core/%.h
		@- mkdir -p obj
		$(CC) $(CFLAGS) $(INCS) -c -o $@ $<

$(PE_OBJS):	obj/%.o : src/ParsEval/%.c inc/ParsEval/%.h
		@- mkdir -p obj
		$(CC) $(CFLAGS) $(INCS) -c -o $@ $<

$(VN_OBJS):	obj/%.o : src/VAnG/%.c inc/VAnG/%.h
		@- mkdir -p obj
		$(CC) $(CFLAGS) $(INCS) -c -o $@ $<
		
$(PE_EXE):	src/ParsEval/parseval.c $(AGN_OBJS) $(PE_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) -fopenmp $(INCS) -o $@ $(AGN_OBJS) $(PE_OBJS) src/ParsEval/parseval.c $(LDFLAGS)

$(CN_EXE):	src/canon-gff3.c $(AGN_OBJS) obj/PeNodeVisitor.o
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/canon-gff3.c obj/PeNodeVisitor.o $(LDFLAGS)
		
$(VN_EXE):	src/VAnG/vang.c $(VN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(VN_OBJS) src/VAnG/vang.c $(LDFLAGS)
