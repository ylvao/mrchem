import numpy as np
import argparse
import basis_set_exchange as bse
import mendeleev
from mendeleev.econf import ElectronicConfiguration

from pyscf import gto, scf


def uncontract_dfunctions(basname='3-21g', atomic_symbol='He'):
    """Fetch 3-21G from basis set exchange and uncontract d-functions if present."""
    # Fetch basis set in NWChem format from BSE.
    bas = bse.get_basis(name=basname, elements=[atomic_symbol], fmt='nwchem', header=True)
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


def get_basis_in_mrchem_format(atomic_symbol):
    bas = bse.convert_formatted_basis_str(uncontract_dfunctions(atomic_symbol=atomic_symbol), 'nwchem', 'dalton')
    atom = mendeleev.element(atomic_symbol)
        
    Z = atom.atomic_number
    if Z > 20:
        nshells = 3
    elif Z > 2:
        nshells = 2
    elif Z > 0:
        nshells = 1
    else:
        raise ValueError('Invalid atomic number!')
    nfuncs = ' '.join(['1' for _ in range(nshells)])
        

    # Ensure correct formatting for MRChem
    new = []
    for line in bas.splitlines():
        if line.strip().startswith('!'):
            continue
        elif line.strip() == '':
            continue
        elif line.strip().startswith('a'):
            continue
        elif line.strip().startswith('H'):
            new.append(' '.join(line.split()[1:]))
        else:
            new.append(' '.join(list(map(lambda x: f'{float(x):.12e}', line.split()))))

    b = []
    b.append('Gaussian basis 3-21G\n')
    b.append('1\n')
    b.append(f'{Z}. 1 {nshells} {nfuncs}\n')
    b.append(f'{atomic_symbol} {0:.12f} {0:.12f} {0:.12f}\n')
    for line in new:
        b.append(line + '\n')

    return "".join(b)


def sad(atomic_symbol):
    """Compute density matrix for passed atom, and write density and basis set to file."""
    # Build molecule input
    basis = uncontract_dfunctions(atomic_symbol=atomic_symbol, basname='3-21g')
    mol = gto.Mole(atom=f'{atomic_symbol} 0.0 0.0 0.0', basis=basis, symmetry=False, verbose=False)
    mol.spin = ElectronicConfiguration(mendeleev.element(atomic_symbol).econf).unpaired_electrons()
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

    # Flatten density matrix and get basis set in mrchem format
    D_fmt = [f'{x:.12e}' for _ in D.tolist() for x in _]
    basis_mrc = get_basis_in_mrchem_format(atomic_symbol)

    # Write files
    with open(atomic_symbol + '.dens', 'w') as d, open(atomic_symbol + '.bas', 'w') as b:
        d.write(f'{dim}\n')
        for line in D_fmt:
            d.write(line + '\n')

        b.write(basis_mrc)
        
    return energy
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate SAD guess files for MRChem')
    parser.add_argument('element', help='Atomic symbol for which you want to generate SAD files.')
    args = parser.parse_args()
    
    e = sad(args.element)
    print(f'{"UHF/3-21G for":15}: {args.element}')
    print(f'{"Total energy":15}: {e:.8f} a.u.')
    print(f'{"Density file":15}: {args.element}.dens')
    print(f'{"Basis set file":15}: {args.element}.bas')