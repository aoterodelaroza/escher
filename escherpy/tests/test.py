#!/bin/bash

nosetests -v --with-coverage --cover-html --cover-package=escherpy

iceweasel cover/index.html
