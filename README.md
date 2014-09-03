Simulation and data analysis code for parametric spiking models of the
cerebellar granule cell layer. Cell and synaptic models necessary to
run the simulations are mantained separately at
https://github.com/epiasini/BillingsEtAl2014_GCL_Models.

Requirements
------------

- Neuron >= 7.3
- neuroConstruct >= 1.7.1
- jNeuroML >= 5.0.1

Python packages:
- `numpy` >= 1.7
- `h5py` >= 2.2
- `networkX` >= 1.9
- `decorator` >= 3.4 (required for networkX)

Note that `networkX` needs to be available both to CPython and to the
Jython interpreter embedded in neuroConstruct. A convenient way of
ensuring this is by installing the module somewhere and by
appropriately modifying the JYTHONPATH environment variable. Same
thing holds for `decorator`.
