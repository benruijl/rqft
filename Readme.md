RQFT
====

RQFT is a program that can compute the renormalization part (i.e., the local counterterm) of *any* Feynman graph up to three loops. It expresses all (sub)counterterms of the graph as massive vacuum bubbles and has built-in IBP reduction tables for these topologies.

Prerequisites
-------------
Install:
-  the latest version of [FORM](https://github.com/vermaseren/form) from the `master` branch
-  [qgraf 3](http://cfif.ist.utl.pt/~paulo/qgraf.html) or higher
-  Python 3.8 or higher


Usage
-----
Specify the process in qgraf using `model.dat` and `qgraf.dat`. Call the output `$process_name.py` and move it one folder up. For this example, we will assume `$process_name` is `ghost_nnnlo`.

Call:
```bash
python generate_process.py ghost_nnnlo
```

This creates a folder `ghost_nnnlo` that contains a file `ghost_nnnlo_$i_in.h` for every graph `$i`, describing its ultraviolet structure. In the folder, running

```bash
python run_process.py
```
will produce the renormalization parts for all graphs in the process. They can be found in the files `ghost_nnnlo_$i_out.h`, where `$i` is the graph index.

Using `Makefile` in that folder generates visualisations of the graph using Mathematica.

To verify 

Extensions
--------------------------
To extend this code to higher loops, the only ingredient that is needed is the reduction table and master evaluations for single-scale massive vacuum bubbles. These can be added to `integrateduv.frm`.

By default, the QCD Feynman rules are provided. To go beyond QCD, the extra Feynman rules can be implemented in `feynman_rules.h`.