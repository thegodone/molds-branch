from rdkit import Chem
from rdkit.Chem import AllChem


def outputmolds(mol,operator,filename):
  numAtoms = mol.GetConformer().GetNumAtoms()
  with open(filename, 'w') as f:
      f.write('THEORY\n'+operator+'\nTHEORY_END\n\n')
      f.write('SCF\n max_iter 50\n rms_density 0.000001\n damping_thresh 1.0\n damping_weight 0.0\n diis_num_error_vect 5\n diis_start_error 0.1\n diis_end_error 0.00000002\nSCF_END\n\n')
      f.write('OPTIMIZATION\n   method gediis\n   total_steps 50\n   electronic_state 0\n   max_gradient 0.00045\n   rms_gradient 0.00030\nOPTIMIZATION_END\n\n')
      f.write('GEOMETRY\n')
      for i in xrange(numAtoms):
          f.write("{:2s}  {:7.4f} {:7.4f} {:7.4f}\n".format(
			  mol.GetAtomWithIdx(i).GetSymbol(),
              mol.GetConformer().GetAtomPosition(i)[0],
              mol.GetConformer().GetAtomPosition(i)[1],
              mol.GetConformer().GetAtomPosition(i)[2]))
      f.write('GEOMETRY_END\n')
      f.close()
  return


# main script to convert the smile to molds file
m = Chem.MolFromSmiles('Cc1ccccc1')
m = Chem.AddHs(m)
AllChem.EmbedMolecule(m,AllChem.ETKDG())
# print the molds format from a smile
outputmolds(m,'am1','toluene_AM1_opt_gediis.in')
