top_srcdir = ../..
srcdir     = .

SHELL      = /bin/sh

# Installation targets
#
prefix      = /home/RiboVision3/R2DT1_4/infernal/infernal-1.1.5/infernal
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

MANS =  infernal\
	cmalign\
	cmbuild\
	cmcalibrate\
	cmconvert\
	cmemit\
	cmfetch\
	cmpress\
	cmsearch\
	cmscan\
	cmstat
# don't install these, since the programs aren't installed
#	sse_cmsearch\
#	cmsearch_dual

.PHONY:  install uninstall clean distclean

install:
	for file in ${MANS}; do \
	   ${INSTALL} -m 0664 $$file.man ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

uninstall:
	for file in ${MANS}; do \
	   rm ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

clean:
	-rm -f *~ 
ifndef V
	@echo '     ' CLEAN man
endif

distclean:
	-rm -f *~ 
	-rm -f *.man
	-rm -f Makefile
