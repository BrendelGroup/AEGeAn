#---Begin configuration---#
#---Begin configuration---#

prefix=/usr/local
GT_INSTALL_DIR=$(prefix)
AGN_LINK=https://github.com/standage/AEGeAn
AGN_DATE=2013
AGN_VERSION=0.9.4-rc

#----End configuration----#
#----End configuration----#

# Binaries
PE_EXE=bin/parseval
CN_EXE=bin/canon-gff3
VN_EXE=bin/vang
LP_EXE=bin/locuspocus
UT_EXE=bin/unittests
BINS=$(PE_EXE) $(CN_EXE) $(VN_EXE) $(LP_EXE)

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
CFLAGS=-Wall -DAGN_DATA_PATH='"$(prefix)/share/aegean"' -fopenmp # -Wno-unused-result
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
INCS=-I $(GT_INSTALL_DIR)/include/genometools/ -I inc/core -I inc/ParsEval -I inc/VAnG -I /usr/include/cairo/ -I /sw/include/cairo/

# Targets
all:		$(BINS) $(UT_EXE) libaegean.a
		

install:	all
		@- test -d $(prefix)/bin || mkdir $(prefix)/bin
		cp $(BINS) $(prefix)/bin/.
		cp libaegean.a $(prefix)/lib/.
		@- test -d $(prefix)/include/aegean || mkdir $(prefix)/include/aegean
		@- rm -f $(prefix)/include/aegean/*
		cp inc/core/*.h $(prefix)/include/aegean/.
		@- test -d $(prefix)/share || mkdir $(prefix)/share
		@- test -d $(prefix)/share/aegean || mkdir $(prefix)/share/aegean
		cp -r data/share/* $(prefix)/share/aegean/.

uninstall:	
		rm -r $(prefix)/$(PE_EXE)
		rm -r $(prefix)/share/parseval

clean:		
		rm -f $(BINS) $(UT_EXE) libaegean.a $(CLSS_MDL_OBJS) inc/core/AgnVersion.h

$(AGN_OBJS):	obj/%.o : src/core/%.c inc/core/%.h inc/core/AgnVersion.h
		@- mkdir -p obj
		$(CC) $(CFLAGS) $(INCS) -c -o $@ $<

$(PE_OBJS):	obj/%.o : src/ParsEval/%.c inc/ParsEval/%.h inc/core/AgnVersion.h
		@- mkdir -p obj
		$(CC) $(CFLAGS) $(INCS) -c -o $@ $<

$(VN_OBJS):	obj/%.o : src/VAnG/%.c inc/VAnG/%.h inc/core/AgnVersion.h
		@- mkdir -p obj
		$(CC) $(CFLAGS) $(INCS) -c -o $@ $<
		
$(PE_EXE):	src/ParsEval/parseval.c $(AGN_OBJS) $(PE_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) $(PE_OBJS) src/ParsEval/parseval.c $(LDFLAGS)

$(CN_EXE):	src/canon-gff3.c $(AGN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/canon-gff3.c $(LDFLAGS)
		
$(VN_EXE):	src/VAnG/vang.c $(VN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(VN_OBJS) src/VAnG/vang.c $(LDFLAGS)

$(LP_EXE):	src/locuspocus.c $(AGN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/locuspocus.c $(LDFLAGS)

$(UT_EXE):	test/unittests.c $(AGN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) test/unittests.c $(LDFLAGS)

libaegean.a:	$(AGN_OBJS)
		ar ru libaegean.a $(AGN_OBJS)

inc/core/AgnVersion.h:	
			@- bash -c "if [ -d .git ]; then perl data/scripts/version.pl > inc/core/AgnVersion.h; else perl data/scripts/version.pl --link=$(AGN_LINK) --date=$(AGN_DATE) --version=$(AGN_VERSION) > inc/core/AgnVersion.h; fi"

test:		$(BINS) $(UT_EXE)
		@- bin/unittests
		@- echo AEGeAn Functional Tests
		@- test/AT1G05320.sh
		@- test/FBgn0035002.sh
