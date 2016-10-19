This is the FDFIT package, for working with DNA force-extension data in MATLAB.


GETTING STARTED
===============

In the "Doc" directory, you'll find a number of sample scripts. Start with
"Tutorial" for the basics. A formatted version of the sample scripts (including
figures) can be found in the "Doc\html" directory ("Tutorial.html", etc.)


REQUIREMENTS
============

MATLAB R2012a or higher (tested with R2014b), including the following Toolboxes:

 - Curve Fitting Toolbox;
 - Optimization Toolbox;
 - Parallel Computing Toolbox (optional).


CREDITS
=======

This code was written by Onno Broekmans, during his PhD research in the lab of
prof. Gijs J.L. Wuite [1]. Prof. Wuite is the corresponding author on the
publication (see "CITATIONS.txt").

The FDFIT package includes some additional libraries, under the "Lib" directory:

 - cursors.m, by T. Montagnon;
 - GUI Layout Toolbox, by Ben Tordoff / The Mathworks;
 - TDMSReader, by Jim Hokanson.

See the respective README/LICENSE files for details.

Thanks to these people for making their code available for re-use!


PERFORMANCE NOTES
=================

The analyses in this package have been run on a pretty powerful workstation
(2x 8-core Intel Xeon CPU, 32 GB RAM), so not a lot of effort has been invested
in optimizing the code for performance. Reproducing the figures from the paper
is therefore likely to require at least an overnight run on an 'average' PC.

The code does make extensive use of parallel processing (MATLAB's "parfor"
loop). You will need the "Parallel Computing Toolbox" for this to work properly
(highly recommended). Without that Toolbox, the "parfor" loops will fall back to
a normal "for" loop, which is much slower.


REPRODUCING THE FIGURES FROM THE PAPER
======================================

Reproducing all the figures from the paper's main text and Supplementary
Information should be as easy as opening this folder in MATLAB, and running:

    >> init
    >> reproducePaper

Then follow the on-screen instructions.

NOTE: You will need to have downloaded the paper's raw data from Figshare [2]
for this to work.


REFERENCES
==========

 [1]: http://www.nat.vu.nl/~gwuite/

 [2]: O.D. Broekmans, G.A. King, G.J. Stephens, and G.J.L. Wuite,
      "Double-stranded DNA force-extension data, as a function of magnesium 
      concentration." http://dx.doi.org/10.6084/m9.figshare.1299526 (2015).

