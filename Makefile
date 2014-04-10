#---Begin configuration---#
#---Begin configuration---#

prefix=/usr/local
AGN_LINK=https://github.com/standage/AEGeAn
AGN_DATE=2014
AGN_VERSION=0.9.4-rc

#----End configuration----#
#----End configuration----#

# Binaries
PE_EXE=bin/parseval
CN_EXE=bin/canon-gff3
VN_EXE=bin/vang
LP_EXE=bin/locuspocus
XT_EXE=bin/xtractore
RP_EXE=bin/pmrna
UT_EXE=bin/unittests
BINS=$(PE_EXE) $(CN_EXE) $(VN_EXE) $(LP_EXE) $(XT_EXE) $(RP_EXE)

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
        -Lsrc/genometools/lib
ifdef lib
  LDFLAGS += -L$(lib)
endif
INCS=-I src/genometools/src -I src/genometools/obj -I /usr/local/include/genometools \
     -I inc/core -I inc/ParsEval -I inc/VAnG       \
     -I /usr/include/cairo -I /sw/include/cairo
ifeq ($(memcheck),yes)
  MEMCHECK=valgrind --leak-check=full --show-reachable=yes --error-exitcode=1 \
                    --suppressions=data/share/libpixman.supp
  MEMCHECKFT=memcheck
endif
LDPATH=LD_LIBRARY_PATH=src/genometools/lib DYLD_LIBRARY_PATH=src/genometools/lib

# Targets
all:		gt agn
		

agn:		$(LP_EXE) $(XT_EXE) $(UT_EXE) libaegean.a
		

install:	all gt-install
		

agn-install:	agn
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

all-clean:	clean gt-clean
		

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

$(XT_EXE):	src/xtractore.c $(AGN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/xtractore.c $(LDFLAGS)

$(RP_EXE):	src/pmrna.c $(AGN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) src/pmrna.c $(LDFLAGS)

$(UT_EXE):	test/unittests.c $(AGN_OBJS)
		@- mkdir -p bin
		$(CC) $(CFLAGS) $(INCS) -o $@ $(AGN_OBJS) test/unittests.c $(LDFLAGS)

libaegean.a:	$(AGN_OBJS)
		ar ru libaegean.a $(AGN_OBJS)

inc/core/AgnVersion.h:	
			@- bash -c "if [ -d .git ]; then perl data/scripts/version.pl > inc/core/AgnVersion.h; else perl data/scripts/version.pl --link=$(AGN_LINK) --date=$(AGN_DATE) --version=$(AGN_VERSION) > inc/core/AgnVersion.h; fi"

test:		all agn-test
		

agn-test:	agn
		@- $(LDPATH) $(MEMCHECK) bin/unittests
		@- $(LDPATH) $(MEMCHECK) bin/locuspocus --outfile=/dev/null data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@- $(LDPATH) $(MEMCHECK) bin/locuspocus --outfile=/dev/null --intloci data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@- $(LDPATH) $(MEMCHECK) bin/locuspocus --outfile=/dev/null --intloci --skipends data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@- $(LDPATH) $(MEMCHECK) bin/locuspocus --outfile=/dev/null --intloci --skipends --verbose data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@- $(LDPATH) $(MEMCHECK) bin/locuspocus --outfile=/dev/null --intloci --skipends --verbose --idformat=GrapeLocus%03lu data/gff3/grape-refr.gff3 data/gff3/grape-pred.gff3
		@- echo AEGeAn Functional Tests
		@- test/xtractore-ft.sh $(MEMCHECKFT)
		@- #test/AT1G05320.sh
		@- #test/FBgn0035002.sh
		@- #test/iLocusParsing.sh

gt:		
		cd src/genometools; make $(GTFLAGS)

gt-install:	
		cd src/genometools; make $(GTFLAGS) install

gt-clean:	
		cd src/genometools; make clean
