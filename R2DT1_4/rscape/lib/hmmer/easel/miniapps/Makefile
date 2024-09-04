top_srcdir = ../..
srcdir     = .

SHELL      = /bin/sh

# Package information
#
PACKAGE         = HMMER
PACKAGE_VERSION = 3.3.2
PACKAGE_TARNAME = hmmer
BASENAME        = hmmer
RELEASE         = 0.48
RELEASEDATE     = "Nov 2020"
COPYRIGHT       = "Copyright (C) 2020 Howard Hughes Medical Institute."
LICENSELINE     = "Freely distributed under the BSD open source license."

## Installation targets
#
prefix      = /home/RiboVision3/R2DT1_4/rscape
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin
libdir      = ${exec_prefix}/lib
includedir  = ${prefix}/include
mandir      = ${datarootdir}/man
docdir      = ${datarootdir}/doc/${PACKAGE_TARNAME}
pdfdir      = ${docdir}
mandir      = ${datarootdir}/man
man1dir     = ${mandir}/man1
man1ext     = .1

INSTALL     = /bin/install -c

## your compiler and compiler flags.
#
CC             = gcc
CFLAGS         = -O3
PTHREAD_CFLAGS = -pthread 
CPPFLAGS       = 
LDFLAGS        = 
DEFS           = -DHAVE_CONFIG_H
LIBS           = -leasel   -lpthread -lm

## list of the miniapps to compile.
#
PROGS = esl-afetch\
	esl-alimanip\
	esl-alimap\
	esl-alimask\
	esl-alimerge\
	esl-alipid\
	esl-alirev\
	esl-alistat\
	esl-compalign\
	esl-compstruct\
	esl-construct\
	esl-histplot\
	esl-mask\
	esl-mixdchlet\
	esl-reformat\
	esl-selectn\
	esl-seqrange\
	esl-seqstat\
	esl-sfetch\
	esl-shuffle\
	esl-ssdraw\
        esl-translate\
	esl-weight

SUBCMDOBJS = \
	cmd_alistat.o    \
	cmd_downsample.o \
	cmd_filter.o     \
	cmd_index.o

# beautification magic stolen from git 
#
QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR  $$subdir; \
		        ${MAKE} -s -C $$subdir
endif


.PHONY: all dev tests check install uninstall clean distclean 

all: 	${PROGS} easel
dev:    ${PROGS} easel
tests:  ${PROGS} easel
check:  ${PROGS} easel

easel:  % : %.c ../libeasel.a ${SUBCMDOBJS}
	${QUIET_GEN}${CC} ${CPPFLAGS} ${CFLAGS} ${PTHREAD_CFLAGS} ${DEFS} ${LDFLAGS} -L.. -I. -I.. -I${srcdir} -I${srcdir}/.. -o $@ $< ${SUBCMDOBJS} ${LIBS}

${PROGS}: % : %.c ../libeasel.a 
	${QUIET_GEN}${CC} ${CPPFLAGS} ${CFLAGS} ${PTHREAD_CFLAGS} ${DEFS} ${LDFLAGS} -L.. -I. -I.. -I${srcdir} -I${srcdir}/.. -o $@ $< ${LIBS}

${SUBCMDOBJS}: %.o : %.c ../libeasel.a
	${QUIET_CC}${CC} -I. -I.. -I${srcdir} -I${srcdir}/.. ${CPPFLAGS} ${CFLAGS} ${PTHREAD_CFLAGS} ${SIMD_CFLAGS} ${DEFS} -c $<

install: ${PROGS} 
	${INSTALL} -m 0755 easel ${DESTDIR}${bindir}/
	for file in ${PROGS}; do \
	   ${INSTALL} -m 0755 $$file ${DESTDIR}${bindir}/ ;\
	done
	for file in ${PROGS}; do \
	   ${INSTALL} -m 0644 $$file.man ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

uninstall:
	rm ${DESTDIR}${bindir}/easel
	for file in $(PROGS); do \
	   rm ${DESTDIR}${bindir}/$$file ;\
	done
	for file in ${PROGS}; do \
	   rm ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

clean:
	-rm -f ${PROGS} easel
	-rm -f *.o *~
	-rm -f *.gcno *.gcda *.gcov
	-rm -f TAGS
	for prog in ${PROGS}; do\
	   if test -d $$prog.dSYM; then rm -rf $$prog.dSYM; fi;\
	done
ifndef V
	@echo '     ' CLEAN miniapps
endif


distclean: clean
	-rm -f ${PROGS}
	-rm -f *.man
	-rm -f Makefile


