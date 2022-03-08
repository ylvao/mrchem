import numpy as np
import basis_set_exchange as bse
import mendeleev
from mendeleev.econf import ElectronicConfiguration

from pyscf import gto, scf


def uncontract_dfunctions(basname='3-21g', atomic_symbol='He'):
    """Fetch 3-21G from basis set exchange and uncontract d-functions if present."""
    # Fetch basis set in NWChem format from BSE.
    bas = bse.get_basis(name=basname, elements=[atomic_symbol], fmt='nwchem', header=False)
    atom = mendeleev.element(atomic_symbol)
    
    # No need to uncontract d-functions when there are none
    # Note: For polarized basis sets this will not work
    # (3-21G does not contain polarization functions)
    if atom.atomic_number < 21:
        return bas
    
    # Get the d function exponents
    # (not very pretty)
    lines = bas.splitlines()
    d_exp = []
    for i, line in enumerate(lines):
        if line.strip().startswith(atomic_symbol) and line.split()[1] == 'D':
            for j, subline in enumerate(lines[i:]):
                if atomic_symbol in subline:
                    continue
                elif 'END' in subline:
                    break
                else:
                    d_exp.append(float(subline.split()[0]))
            break
    
    # Build the custom basis set
    new = [line for line in lines if 'D' not in line.split() and 'END' not in line][:-len(d_exp)]
    for e in d_exp:
        new.append(f'{atomic_symbol}    D')
        new.append(f'{" "*(len(atomic_symbol)+4)}{e:.10f}       {1:.10f}')
    new.append('END')
        
    return '\n'.join(new)


def sad(atom):
    """Compute density matrix for passed atom, and write density and basis set to file."""
    # Build molecule input
    basis = uncontract_dfunctions(atomic_symbol=atom, basname='3-21g')
    mol = gto.Mole(atom=f'{atom} 0.0 0.0 0.0', basis=basis, symmetry=False, verbose=False)
    mol.spin = ElectronicConfiguration(mendeleev.element(atom).econf).unpaired_electrons()
    mol.build()
    
    # Set SCF method
    qc = scf.UHF(mol)
    
    # Set some SCF parameters
    qc.conv_tol = 1e-10
    qc.max_cycle = 500
    
    # Perform the calculation
    energy = qc.kernel()
    
    # Compute AO density matrix with some thresholding
    Da, Db = qc.make_rdm1()
    D = Da + Db
    zeros = np.abs(D) <= 1e-13
    D[zeros] = 0.0
    dim = D.shape[0]

    # Flatten density matrix and get basis set in dalton format
    D_fmt = [f'{x:.12e}' for _ in D.tolist() for x in _]
    dalton_basis = bse.convert_formatted_basis_str(basis, 'nwchem', 'dalton')

    # Write files
    with open(atom + '.dens', 'w') as d, open(atom + '.bas', 'w') as b:
        d.write(f'{dim}\n')
        for line in D_fmt:
            d.write(line + '\n')

        b.write(dalton_basis)
        
    return energy
        
        
if __name__ == '__main__':
    row5 = [el.symbol for el in mendeleev.get_all_elements() if el.period == 5]
    for element in row5:
        e = sad(element)
        print(f'{element:5}: {e:.8f} Hartree')
