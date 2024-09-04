#
# Kinfold subpackage
#

AC_DEFUN([RNA_ENABLE_PKG_KINFOLD],[
  RNA_ADD_PACKAGE([kinfold],
                  [Kinfold program],
                  [yes],[],[],
                  [${srcdir}/src/Kinfold/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([kinfold],[
    AC_CONFIG_SUBDIRS([src/Kinfold])
  ])

  AM_CONDITIONAL(MAKE_KINFOLD, test "x$with_kinfold" != "xno")
])


#
# RNAforester subpackage
#
AC_DEFUN([RNA_ENABLE_PKG_FORESTER],[
  RNA_ADD_PACKAGE([forester],
                  [RNAforester program],
                  [yes],[],[],
                  [${srcdir}/src/RNAforester/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([forester],[
    AC_CONFIG_SUBDIRS([src/RNAforester])
  ])

  AM_CONDITIONAL(MAKE_FORESTER, test "x$with_forester" != "xno")
])


#
#
#

AC_DEFUN([RNA_ENABLE_PKG_CLUSTER],[
  RNA_ADD_PACKAGE([cluster],
                  [AnalyseSeqs and AnalyseDists],
                  [no],[],[],
                  [${srcdir}/src/Cluster/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([cluster],[
    AC_DEFINE([WITH_CLUSTER], [1], [Analyse{Dists,Seqs}])
    AC_SUBST([CLUSTER_DIR], [Cluster])
    AC_CONFIG_FILES([src/Cluster/Makefile])
  ])

  AM_CONDITIONAL(MAKE_CLUSTER, test "x$with_cluster" = "xyes")
])

#
# Kinwalker subpackage
#

AC_DEFUN([RNA_ENABLE_PKG_KINWALKER],[
  RNA_ADD_PACKAGE([kinwalker],
                  [Kinwalker program],
                  [no],[],[],
                  [${srcdir}/src/Kinwalker/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([kinwalker],[
    AC_CONFIG_SUBDIRS([src/Kinwalker])
  ])

  AM_CONDITIONAL(MAKE_KINFOLD, test "x$with_kinwalker" != "xyes")
])


#
# RNAlocmin subpackage
#
AC_DEFUN([RNA_ENABLE_PKG_RNALOCMIN],[
  RNA_ADD_PACKAGE([rnalocmin],
                  [RNAlocmin program],
                  [yes],[],[],
                  [${srcdir}/src/RNAlocmin/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([rnalocmin],[
    AC_CONFIG_SUBDIRS([src/RNAlocmin])
  ])

  AM_CONDITIONAL(MAKE_RNALOCMIN, test "x$with_rnalocmin" != "xno")
])


