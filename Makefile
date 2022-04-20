BLK         =   black
NOSE        =   nose2
PIP         =   pip

BASEDIR     =   $(CURDIR)
SRCDIR      =   ${BASEDIR}/skymap

help:
	@echo 'Makefile for skymap'
	@echo 'Usage:'
	@echo 'make black           reformat the code using black code formatter'
	@echo 'make clean           remove temporary files'
	@echo 'make install         install the package locally'
	@echo 'make tests           run the unit tests'

black:
	${BLK} *.py */*.py

clean:
	rm -f ${SRCDIR}/*.pyc
	rm -rf ${SRCDIR}/__pycache__
	rm -rf ${BASEDIR}/*/__pycache__
	rm -rf ${BASEDIR}/build
	rm -rf ${BASEDIR}/dist
	rm -rf ${BASEDIR}/skymap.egg-info

install:
	${PIP} install .

tests:
	${NOSE}

.PHONY: help black clean install tests