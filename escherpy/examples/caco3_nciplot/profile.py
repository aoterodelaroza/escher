#!/usr/bin/env python
# coding: utf-8


from memory_profiler import profile


@profile
def main():
    import escherpy as esc
    mol = esc.Molecule()

    mol.structfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
    mol.readstruct()
    mol.stickball()

    densfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-dens.cube'
    gradfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
    mol.nciplot(densfile, gradfile)
    mol.nciplot2D(densfile, gradfile)

    mol.show()

main()

