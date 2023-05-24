# GraphCutOptimization

from https://vision.cs.uwaterloo.ca/code/ section "Multi-label optimization" :

<div style="flex">
    <img src="Tsukuba_tiny.png">
    <img src="Baseball_tiny.png">
    <img src="Oxford_tiny.jpg">
</div>

The **gco-v3.0** library is for optimizing multi-label energies via the α-expansion and α-β-swap algorithms. It supports energies with any combination of unary, pairwise, and label cost terms. Written in C++, ~~it comes bundled with a MATLAB wrapper~~. The code was developed by [Olga Veksler](http://mouse.cs.uwaterloo.ca/viki/Olga_Veksler) and [Andrew Delong](http://mouse.cs.uwaterloo.ca/viki/Andrew_Delong). Andreas Mueller (U. Bonn) has written a [Python wrapper for gco](http://peekaboo-vision.blogspot.ca/2012/05/graphcuts-for-python-pygco.html).

**New in version 3.0** – older versions [here](http://www.csd.uwo.ca/~olga/OldCode.html)

- *Label cost* support in α-expansion, including costs on subsets of labels ("category costs").
- *Sparse data cost* support, for when each label is feasible for only a small fraction of sites; a fast special case!
- *Adaptive cycles* for α-expansion are used by default; often faster, but sometimes slower, so try both.

**Downloads:**

- Source code: [gco-v3.0.zip](https://vision.cs.uwaterloo.ca/files/gco-v3.0.zip) – Apr 28 2010 – patched Oct 14, 2014 to allow float/double energy terms
- Relevant papers: [PAMI2001](http://www.csd.uwo.ca/faculty/olga/Papers/pami01_final.pdf), [PAMI2004a](http://www.csd.uwo.ca/~yuri/Papers/pami04.pdf), [PAMI2004b](http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/papers/KZ-PAMI-graph_cuts.pdf), and (if using label costs) [IJCV2012](http://www.csd.uwo.ca/~yuri/Papers/ijcv10_labelcost.pdf).

**Requirements:**

- Visual C++ 2005 (VC8); GCC 4.03 (warning: not well-tested with GCC)

- ~~MATLAB 7.4 (R2007a) for 32-bit wrapper; MATLAB 7.6 (R2008) for 64-bit wrapper~~ (I removed the MATLAB code)

---

Original README is `GCO_README.TXT`