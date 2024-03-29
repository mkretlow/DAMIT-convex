# DAMIT-convex
Asteroid light curve inversion code by Kaasalainen and Durech from DAMIT [website](https://astro.troja.mff.cuni.cz/projects/damit/pages/software_download).

## Description
The source codes of light-curve inversion routines together with brief manuals, example lightc urves, and the code for the direct problem are available. The code was developed by *Mikko Kaasalainen* in Fortran and converted to C by *Josef Durech*. There are two programs for light-curve inversion: **convexinv**, that optimizes all parameters and uses spherical harmonics functions for shape representation, and **conjgradinv**, that optimizes only shape and uses directly facet areas as parameters – it should be used at the final stage of the inversion process for 'polishing' the final shape model. You can also compute synthetic light curves with lcgenerator.

## License
This software is licensed under [CC Attribution 4.0 international License](https://creativecommons.org/licenses/by/4.0/legalcode).
