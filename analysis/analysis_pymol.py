# !python
# Visualize DiffDock results

import pathlib,os
from pymol import cmd,util
import argparse
from matplotlib.cm import get_cmap

def parse_args():
  parser = argparse.ArgumentParser(description='Visualize DiffDock results with PyMOL')
  parser.add_argument('--protein_pdb', action="store", help='Input Target Structure',type=str, required=True)
  parser.add_argument('--ligand_dir', action="store", help='Directory of DiffDock Prediction',type=str, required=True)
  parser.add_argument('--ligand_name', action="store", help='Input Ligand Name',type=str, required=True)
  parser.add_argument('--rank', action="store",  nargs=2,default= [1, 20], help="Show number of rank")
  parser.add_argument('--suffix', action="store", default='diffdock', help="Show number of rank")
  parser.add_argument('--cmap', action="store", default='bwr_r',type=str, help="Colormap for ligand")
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

  cmap = get_cmap(args.cmap)
  def ligand_rank_color(cmap,data_list):
    num_color=cmap.N
    data_min,data_max=min(data_list),max(data_list)

    # reflect to a scaled data list
    colors=[list(cmap(int(num_color * ((item-data_min)/(data_max-data_min)))))[:3] for item in data_list]

    return colors



  # confidences

  sdf_confidence=[x for x in os.listdir(f'{LIGAND_FP}') if '_confidence' in x and x.endswith('.sdf') and x.startswith('rank')]
  num_predicts = len(sdf_confidence)
  print(f'{num_predicts}: {sdf_confidence}')

  # parse confidence score from sdf filename
  confidence_dict={}
  for rank in range(1,num_predicts+1):
    sdf_rank=[x for x in sdf_confidence if x.startswith(f'rank{rank}_confidence')][0]
    confidence_score=float(sdf_rank.split('confidence')[1].split('.sdf')[0])
    print(f'reading confidence {confidence_score} from rank {rank} file: {sdf_rank}')
    confidence_dict[rank]=confidence_score

  # fetch color list
  rank_color = ligand_rank_color(cmap=cmap, data_list=confidence_dict.values())

  cmd.load(filename=PROTEIN_FP,object=PROTEIN_NAME)

  for i in range(int(RANK_NUM[0]),int(RANK_NUM[1])+1):
    if os.path.exists(f'{LIGAND_FP}/rank{i}_reverseprocess.pdb'):
      # cmd.load( filename [,object [,state [,format [,finish [,discrete [,multiplex ]]]]]] )
      cmd.load(filename=f'{LIGAND_FP}/rank{i}_reverseprocess.pdb',
               object=f'{LIGAND_NAME}_{i}',
               state=2,
               #discrete=1,
               multiplex=0,
               )
      cmd.do(f'h_fix {LIGAND_NAME}_{i}')
      print(f'color rank#{i} with color {rank_color[i-1]}')
      cmd.set_color(f'color_ranked_{i}', rank_color[i-1])
      cmd.color(f'color_ranked_{i}', f'({LIGAND_NAME}_{i})')
      util.cnc(f'({LIGAND_NAME}_{i})', _self=cmd)

    else:
      print(f'file not found: {LIGAND_FP}/rank{i}_reverseprocess.pdb')


  cmd.set('surface_cavity_mode',4)
  cmd.set('transparency',0.5)
  cmd.set('surface_color','gray70')
  cmd.show('surface',PROTEIN_NAME)
  cmd.hide('cartoon')


  cmd.set('stick_color','firebrick',PROTEIN_NAME)


  cmd.center(f'{PROTEIN_NAME} and hetatm')

  os.makedirs('pse',exist_ok=True)
  cmd.save(f'pse/{PROTEIN_NAME}_{LIGAND_NAME}_{args.suffix}.pze')