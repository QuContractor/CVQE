from openfermion.hamiltonians import MolecularData
from openfermionpyscf import run_pyscf

from openfermion.transforms import get_fermion_operator, get_sparse_operator, \
jordan_wigner,bravyi_kitaev
import numpy as np




element_names = ['H', 'Li']
basis = 'sto-6g'
charge = 0
multiplicity = 1

# Add points for a full dissociation curve from 0.1 to 3.0 angstroms
n_points=50
spacings = np.linspace(0.5,5,n_points)

# Set run options
#run_scf = 1
#run_mp2 = 1
#run_cisd = 1
#run_ccsd = 1
#run_fci = 1
#verbose = 1


def molecule_data(spacing):
    description = "{}".format(spacing)
    geometry = [[element_names[0], [0, 0, 0]],
                [element_names[1], [0, 0, spacing]]]
    molecule = MolecularData(geometry,
                             basis,
                             multiplicity,
                             charge,
                             description)

    molecule = run_pyscf(molecule)
    return molecule

def partial_average(key,drop_qubits=[1,3,5],avg_dict={'X':0,'Y':0,'Z':1},init_fock_one=[]):
    #  Construct effective Hamiltonian on qubits 0,2,4,
    # by average on |000> of qubit 1,3,5 
    key=list(key)
    factor=1
    new_key=[]
    for k in key:
        if k[0] not in drop_qubits:
            new_key.append(k)
        if k[0] in drop_qubits:
            if k[0] in init_fock_one:
                factor*=-avg_dict[k[1]]
            else:
                factor*=avg_dict[k[1]] # init_fock=1 then \sigma_z=-1
        #print(key)
    new_key=tuple(new_key)
    return (new_key,factor) 

def full_ham_terms(spacing,dead_space = range(1),active_space = range(1,4)):
    molecule=molecule_data(spacing)

    molecular_hamiltonian = molecule.get_molecular_hamiltonian(
    occupied_indices=dead_space,
    active_indices=active_space)

    fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
    qubit_hamiltonian = bravyi_kitaev(fermion_hamiltonian)
    qubit_hamiltonian.compress()
    terms_dict=qubit_hamiltonian.terms
    return terms_dict

def hartree_energy(spacing):
    # Hartree energy   
    terms_dict=full_ham_terms(spacing)

    reduced_terms=[]
    for key in terms_dict.keys():
        rt=partial_average(key,drop_qubits=[0,1,2,3,4,5],init_fock_one=[0])   
        reduced_terms.append(rt) 

    ham_terms=np.array([f[0] for f in reduced_terms])
    factors=np.array([f[1] for f in reduced_terms])
    cs=np.array([c for c in terms_dict.values()])
    cs_rescale=np.multiply(factors,cs)

    reduced_terms_rescale=[]
    for i in range(len(reduced_terms)):
        if cs_rescale[i] !=0:
            reduced_terms_rescale.append((reduced_terms[i][0],cs_rescale[i]))
    reduced_terms_rescale

    sim_dict={}
    for term in reduced_terms_rescale:
        if term[0] not in sim_dict.keys():
            sim_dict[term[0]]=term[1]
        else:
            sim_dict[term[0]]+=term[1]

    hartree_energy=sim_dict[()]
    return hartree_energy

def coeff_ham(spacing):
    terms_dict=full_ham_terms(spacing)
    #reduce 6 qubits to 3qubits.
    reduced_terms=[]
    for key in terms_dict.keys():
        rt=partial_average(key)   
        reduced_terms.append(rt)


    #collect same terms in sim_dict
    ham_terms=np.array([f[0] for f in reduced_terms])
    factors=np.array([f[1] for f in reduced_terms])
    cs=np.array([c for c in terms_dict.values()]) # due to partial average
    cs_rescale=np.multiply(factors,cs)

    reduced_terms_rescale=[]
    for i in range(len(reduced_terms)):
        if cs_rescale[i] !=0:
            reduced_terms_rescale.append((reduced_terms[i][0],cs_rescale[i]))
    
    sim_dict={}
    for term in reduced_terms_rescale:
        if term[0] not in sim_dict.keys():
            sim_dict[term[0]]=term[1]
        else:
            sim_dict[term[0]]+=term[1]

    coeff_h=[v for v in sim_dict.values()]
    ham_term=[v for v in sim_dict.keys()]

    return coeff_h,ham_term
    
if __name__ == "__main__":
    
    coeff_hams=[]
    ham_terms=[]
    hartree_energies=[]
    for i in range(n_points):
        coeff_h,ham_term=coeff_ham(spacings[i])
        coeff_hams.append(coeff_h)
        ham_terms.append(ham_term)
        hartree_energies.append(hartree_energy(spacings[i]))

    coeff_hams=np.array(coeff_hams)
    #for h in ham_terms:
    #    print(h)
    hartree_energies=np.array(hartree_energies)
    np.savetxt("coeff_hams_more.txt",coeff_hams)
    np.savetxt("spacings_more.txt",spacings)
    np.savetxt("hartree_energies_more.txt",hartree_energies)
