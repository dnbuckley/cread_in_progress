#
#  Copyright (C) 2006-2015 Cold Spring Harbor Laboratory
#                          University of Southern California
#                          Authors: Andrew D. Smith, Pavel Sumazin
#                                   and Michael Q. Zhang
#
#  Authors: Andrew D. Smith and Sarah S. Ma
#
#  This file is part of CREAD.
#
#  CREAD is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  CREAD is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with CREAD; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

PROGS = basecomp bbcomp decoymasker getcontext \
	formattable joinmotifs motifclass uniqmotifs \
	matcompare feateval kmercomp sortmotifs \
        storm

CXX = g++
CFLAGS = -Wall -fmessage-length=50
CXXFLAGS = -std=c++11 -L/home/psw-is-psw/bin/gsl/lib -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = ../common
INCLUDEDIRS =  $(SMITHLAB_CPP) $(COMMON_DIR)

INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lgsl -lgslcblas

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin


$(PROGS): $(addprefix $(SMITHLAB_CPP)/, \
	smithlab_os.o smithlab_utils.o GenomicRegion.o OptionParser.o)

bbcomp decoymasker getcontext formattable motifclass kmercomp sortmotifs storm: \
          $(addprefix $(COMMON_DIR)/, FastaFile.o)

getcontext joinmotifs motifclass uniqmotifs matcompare sortmotifs storm: \
          $(addprefix $(COMMON_DIR)/,  Motif.o MotifSite.o Matrix.o Pattern.o)

formattable: $(addprefix $(COMMON_DIR)/, WordTable.o ScoringMatrix.o Matrix.o)

joinmotifs uniqmotifs matcompare: $(addprefix $(COMMON_DIR)/, MatCompMethods.o)

motifclass storm: $(addprefix $(COMMON_DIR)/, SuffixTree.o ScoringMatrix.o)

uniqmotifs matcompare sortmotifs: $(addprefix $(COMMON_DIR)/, \
	Module.o ModuleSite.o)

feateval kmercomp: $(addprefix $(COMMON_DIR)/,  Pattern.o SuffixTree.o )

kmercomp sortmotifs: $(addprefix $(COMMON_DIR)/, Word.o WordMatcher.o WordSite.o)

mafutil: $(addprefix $(COMMON_DIR)/, GenomeUtils.o \
	GenomeAlignment.o BEDFile.o)

sortmotifs: $(addprefix $(COMMON_DIR)/, PatternFactory.o)

storm: $(addprefix $(COMMON_DIR)/, IndexedSuffixTree.o GeneralizedSuffixTree.o WordTable.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
