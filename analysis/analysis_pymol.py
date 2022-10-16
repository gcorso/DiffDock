# !python
# Visualize DiffDock results

import pathlib,os
from pymol import cmd
import argparse


def parse_args():
  parser = argparse.ArgumentParser(description='Visualize DiffDock results with PyMOL')
  parser.add_argument('--protein_pdb', action="store", help='Input Target Structure',type=str, required=True)
  parser.add_argument('--ligand_dir', action="store", help='Directory of DiffDock Prediction',type=str, required=True)
  parser.add_argument('--ligand_name', action="store", help='Input Ligand Name',type=str, required=True)
  parser.add_argument('--rank', action="store",  nargs=2,default= [1, 20], help="Show number of rank")
  parser.add_argument('--suffix', action="store", default='diffdock', help="Show number of rank")
  return parser.parse_args()


if __name__ == '__main__':
  args=parse_args()
  cmd.reinitialize()

  # protein structure w/o ligand for diffdock
  PROTEIN_FP=args.protein_pdb
  LIGAND_FP=args.ligand_dir
  LIGAND_NAME=args.ligand_name

  RANK_NUM=args.rank


  PROTEIN_FP=pathlib.Path(PROTEIN_FP).resolve()
  PROTEIN_NAME=PROTEIN_FP.stem

  LIGAND_FP=pathlib.Path(LIGAND_FP).resolve()

  cmd.load(filename=PROTEIN_FP,object=PROTEIN_NAME)

  for i in range(int(RANK_NUM[0]),int(RANK_NUM[1])+1):
    try:
      # cmd.load( filename [,object [,state [,format [,finish [,discrete [,multiplex ]]]]]] )
      cmd.load(filename=f'{LIGAND_FP}/rank{i}_reverseprocess.pdb',
               state=1,
               object=f'{LIGAND_NAME}_{i}',
               format='pdb',
               discrete=1,
               multiplex=0)
    except:
      print(f'error occurs while loading prediction #{i}.')

  cmd.set('surface_cavity_mode',4)
  cmd.set('transparency',0.5)
  cmd.set('surface_color','gray70')
  cmd.show('surface',PROTEIN_NAME)
  cmd.hide('cartoon')

  cmd.set('stick_color','lime',f'{LIGAND_NAME}_*')
  cmd.set('stick_color','firebrick',PROTEIN_NAME)


  cmd.center(f'{PROTEIN_NAME} and hetatm')

  os.makedirs('pse',exist_ok=True)
  cmd.save(f'pse/{PROTEIN_NAME}_{LIGAND_NAME}_{args.suffix}.pze')