#!/usr/bin/env python
# coding: utf-8

# A multiwfn critical points analyzer


def test_xyz():
    from escherpy import escher_data
    from escherpy.parser import Parser

    parser = Parser()
    parser.xyz(escher_data + 'mol/nitromethane/nitromethane.xyz')

def test_multiwfn():
    from escherpy import escher_data
    from escherpy.parser import Parser

    parser = Parser()
    parser.multiwfn(escher_data + 'mol/nitromethane/CPprop.txt')

def test_cubeheader():
    from escherpy import escher_data
    from escherpy.parser import Parser

    parser = Parser()
    parser.cubeheader(escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube')

def test_cube():
    from escherpy import escher_data
    from escherpy.parser import Parser

    parser = Parser()
    parser.cube(escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube')

def test_basin():
    from escherpy import escher_data
    from escherpy.parser import Parser

    parser = Parser()
    parser.basin(escher_data + 'cryst/aragonite/aragonite-1.basin')
