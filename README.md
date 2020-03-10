# Quantum snake algorithm: 
* This is the source code for the paper: 

Dan-Bo Zhang and Tao Yin, Collective optimization for variational quantum eigensolvers, [Phys. Rev. A **101**, 032311(2020).](https://link.aps.org/doi/10.1103/PhysRevA.101.032311) ([arXiv: 1910.14030](https://arxiv.org/abs/1910.14030))

### Collective variational quantum eigensolvers (CVQE) for quantum chemistry

We propose a hybrid quantum-classical algorithm that can provide collective optimization for the VQE
to solve a group of related Hamiltonians more efficiently.

A brief introduction for VQE can be found in [QuContractor's github repository.](https://github.com/QuContractor/VQE_tutorial/)

### Simulate molecules with varied bond lengths
We simulated different molecules (H2, LiH, HeH+) with CVQE. Numeral simulations show that the CVQE exhibits clear collective behavior in the optimization process of updating parameters. 

* The calculation for H2 can be found in [CVQE_H2 with HiQ/projectQ](https://github.com/QuContractor/CVQE/blob/master/H2/snake_VQE_H2_pjq.ipynb), and [CVQE_H2 with IBM Qiskit](https://github.com/QuContractor/CVQE/blob/master/H2/snake_VQE_H2_IBM.ipynb)
* The calculation for LiH can be found in [CVQE_LiH](https://github.com/QuContractor/CVQE/blob/master/LiH/snake_VQE_LiH.ipynb)
* The calculation for HeH+ can be found in [CVQE_HeH+](https://github.com/QuContractor/CVQE/blob/master/HeH%2B/snake_VQE_HeH.ipynb)

### Avoid local minimum

The snake algorithm gives rise to collective motion of parameters of different tasks
that can avoid being trapped in local minimums. 

The example with H2 can be found in [CVQE_H2](https://github.com/QuContractor/CVQE/blob/master/H2/snake_VQE_H2_2theta.ipynb).



### Cite:

@article{PhysRevA.101.032311,
  title = {Collective optimization for variational quantum eigensolvers},
  author = {Zhang, Dan-Bo and Yin, Tao},
  journal = {Phys. Rev. A},
  volume = {101},
  issue = {3},
  pages = {032311},
  numpages = {8},
  year = {2020},
  month = {Mar},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevA.101.032311},
  url = {https://link.aps.org/doi/10.1103/PhysRevA.101.032311}
}



