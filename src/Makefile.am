SUBDIRS = distributions 

jagsmod_LTLIBRARIES = vonmises.la

vonmises_la_SOURCES = vonmises.cc

vonmises_la_CPPFLAGS = -I$(includedir)/JAGS

vonmises_la_LIBADD = distributions/vonmisesdist.la
			
if WINDOWS
vonmises_la_LIBADD += -ljags-$(JAGS_MAJOR) -ljrmath-0
else
vonmises_la_LIBADD += -ljags -ljrmath
endif

vonmises_la_LDFLAGS = -module -avoid-version
if WINDOWS
vonmises_la_LDFLAGS += -no-undefined
endif
