
.PHONY: help test testcoverage coverage 

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  test            to test all files"
	@echo "  testcoverage    to test all files and visualize tests coverage"
test:
	nosetests -v 

testcoverage:
	nosetests -v --with-coverage --cover-html --cover-package=escherpy
	iceweasel cover/index.html

