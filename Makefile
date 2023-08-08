BLK         =   black
PIP         =   pip3

BASEDIR     =   $(CURDIR)
SRCDIR      =   ${BASEDIR}/skymap

help:
	@echo 'Makefile for skymap'
	@echo 'Usage:'
	@echo 'make black           reformat the code using black code formatter'
	@echo 'make clean           remove temporary files'
	@echo 'make install         install the package locally'
	@echo 'make test            run the non-interactive regression tests'
	@echo 'make testall         run all regression tests'

black:
	${BLK} *.py */*.py */*/*.py

clean:
	rm -f ${SRCDIR}/*.pyc
	rm -rf ${SRCDIR}/__pycache__
	rm -rf ${BASEDIR}/*/__pycache__
	rm -rf ${BASEDIR}/build
	rm -rf ${BASEDIR}/dist
	rm -rf ${BASEDIR}/skymap.egg-info

install:
	${PIP} install .

test:
	pytest --verbose -m 'not interactive'

testall:
	pytest --verbose -s

.PHONY: help black clean install test testall