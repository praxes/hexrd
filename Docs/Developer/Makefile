### Makefile --- for mdef/XRD developer documentation
## Subversion $Id$:

all:  force docs
	cd XRD; make
	cd GUI; make

PYTH_DIR = ../../

include PyDocs.mk

clean: 
	make cleanHTML
	cd XRD; make cleanHTML
	cd GUI; make cleanHTML

force:
