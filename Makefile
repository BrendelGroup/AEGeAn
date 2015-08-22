#---Begin configuration---#
#---Begin configuration---#

prefix=/usr/local

#----End configuration----#
#----End configuration----#

SHELL := bash

# Binaries
PE_EXE=bin/parseval
CN_EXE=bin/canon-gff3
LP_EXE=bin/locuspocus
GV_EXE=bin/gaeval
XT_EXE=bin/xtractore
RP_EXE=bin/pmrna
TD_EXE=bin/tidygff3
UT_EXE=bin/unittests
INSTALL_BINS=$(PE_EXE) $(CN_EXE) $(LP_EXE) $(GV_EXE) $(XT_EXE) $(RP_EXE) $(TD_EXE)
BINS=$(INSTALL_BINS) $(UT_EXE)

#----- Source, header, and object files -----#

# AEGeAn core class and module files
AGN_SRCS=$(wildcard src/core/Agn*.c)
AGN_OBJS=$(patsubst src/core/%.c,obj/%.o,$(AGN_SRCS))
AGN_HDRS=$(patsubst src/core%.c,inc/core/%.h,$(AGN_SRCS))

# Compilation settings
CC=gcc
CFLAGS=-Wall -DAGN_DATA_PATH='"$(prefix)/share/aegean"' -Wno-unused-result
GTFLAGS=prefix=$(prefix)
ifeq ($(cairo),no)
  CFLAGS += -DWITHOUT_CAIRO
  GTFLAGS += cairo=no
endif
ifneq ($(errorcheck),no)
  CFLAGS += -Werror
endif
ifeq ($(optimize),yes)
  CFLAGS += -O3
endif
ifneq ($(64bit),no)
  CFLAGS += -m64
  GTFLAGS += 64bit=yes
endif
ifneq ($(debug),no)
  CFLAGS += -g
endif
LDFLAGS=-lgenometools -lm \
        -L$(prefix)/lib \
        -L/usr/local/lib
ifdef lib
  LDFLAGS += -L$(lib)
endif
INCS=$(shell pkg-config --silence-errors --cflags-only-I cairo) \
     -I inc/core \
     -I /usr/local/include/genometools \
     -I $(prefix)/include/genometools
ifeq ($(memcheck),yes)
  MEMCHECK=valgrind --leak-check=full --show-reachable=yes --error-exitcode=1 \
                    --suppressions=data/misc/libpixman.supp \
                    --suppressions=data/misc/libpango.supp
  MEMCHECKFT=memcheck
endif

# Targets
all:		$(BINS) libaegean.a
		

install:	all
		@ mkdir -p $(prefix)/bin/
		@ mkdir -p $(prefix)/lib/
		@ mkdir -p $(prefix)/include/aegean/
		@ mkdir -p $(prefix)/share/aegean/
		@ rm -f $(prefix)/include/aegean/*
		cp $(INSTALL_BINS) $(prefix)/bin/.
		cp libaegean.a $(prefix)/lib/.
		cp inc/core/*.h $(prefix)/include/aegean/.
		cp -r data/share/* $(prefix)/share/aegean/.

install-scripts:
		@ mkdir -p $(prefix)/bin/
		cp data/scripts/*.p? $(prefix)/bin/.

uninstall:	
		for exe in $(INSTALL_BINS); do rm -r $(prefix)/$$exe; done
		rm -r $(prefix)/include/aegean/
		rm -r $(prefix)/share/aegean/
		rm $(prefix)/lib/libaegean.a

clean:		
		rm -rf $(BINS) libaegean.a $(AGN_OBJS) inc/core/AgnVersion.h bin/*.dSYM

$(AGN_OBJS):	obj/%.o : src/core/%.c inc/core/%.h inc/core/AgnVersion.h
		@- mkdir -p obj
		@ echo "[compile $*]"
		@ $(CC) $(CFLAGS) $(INCS) -c -o $@ $<

$(PE_EXE):	src/ParsEval/parseval.c src/ParsEval/pe_options.c src/ParsEval/pe_utils.c src/ParsEval/pe_options.h src/ParsEval/pe_utils.h $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile ParsEval]"
		@ $(CC) $(CFLAGS) $(INCS) -I src/ParsEval -o $@ $(AGN_OBJS) src/ParsEval/parseval.c src/ParsEval/pe_options.c src/ParsEval/pe_utils.c $(LDFLAGS)

$(CN_EXE):	src/canon-gff3.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile CanonGFF3]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/canon-gff3.c $(LDFLAGS)

$(LP_EXE):	src/locuspocus.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile LocusPocus]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/locuspocus.c $(LDFLAGS)

$(GV_EXE):	src/gaeval.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile GAEVAL]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/gaeval.c $(LDFLAGS)

$(XT_EXE):	src/xtractore.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile Xtractore]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/xtractore.c $(LDFLAGS)

$(RP_EXE):	src/pmrna.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile $@]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/pmrna.c $(LDFLAGS)

$(TD_EXE):	src/tidygff3.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile $@]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/tidygff3.c $(LDFLAGS)

$(UT_EXE):	test/unittests.c $(AGN_OBJS)
		@ mkdir -p bin
		@ echo "[compile unit tests]"
		@ $(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) test/unittests.c $(LDFLAGS)

libaegean.a:	$(AGN_OBJS)
		@ echo "[create libaegean]"
		@ ar ru libaegean.a $(AGN_OBJS)

inc/core/AgnVersion.h:	
			@- echo "[print $@]"
			@ data/scripts/version.py > $@

test:		agn-test
		

agn-test:	all
		@ $(MEMCHECK) bin/unittests
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --skipends data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --skipends --verbose data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --skipends --verbose --idformat=GrapeLocus%03lu data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --verbose data/gff3/dmel-pseudofeat-sort-test-in.gff3 2> >(grep -v 'not unique')
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null data/gff3/amel-plap.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --cds data/gff3/amel-lsm.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --introngenes data/gff3/mrot-cst.gff3
		@ $(MEMCHECK) bin/locuspocus --outfile=/dev/null --introngenes --delta=5 data/gff3/mrot-cst.gff3
		@ $(MEMCHECK) bin/tidygff3 < data/gff3/grape-refr.gff3 > /dev/null
		@ echo AEGeAn Functional Tests
		@ test/AT1G05320.sh $(MEMCHECKFT)
		@ test/FBgn0035002.sh $(MEMCHECKFT)
		@ test/iLocusParsing.sh $(MEMCHECKFT)
		@ test/xtractore-ft.sh $(MEMCHECKFT)
		@ test/canon-gff3-ft.sh $(MEMCHECKFT)
		@ test/gaeval-ft.sh $(MEMCHECKFT)
		@ test/AmelOGSvsNCBI.sh $(MEMCHECKFT)
		@ test/align-convert.sh
		@ test/misc-ft.sh $(MEMCHECKFT)


