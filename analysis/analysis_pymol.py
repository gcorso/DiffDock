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
  protein_fp=args.protein_pdb
  ligand_fp=args.ligand_dir
  ligand_alias=args.ligand_name
  rank_num=args.rank
  protein_fp=pathlib.Path(protein_fp).resolve()
  protein_alias=protein_fp.stem
  ligand_fp=pathlib.Path(ligand_fp).resolve()
  cmap = get_cmap(args.cmap)

  def ligand_rank_color(cmap,data_list):
    num_color=cmap.N
    data_min,data_max=min(data_list),max(data_list)

    # reflect to a scaled data list
    colors=[list(cmap(int(num_color * ((item-data_min)/(data_max-data_min)))))[:3] for item in data_list]
    return colors

  # confidences
  sdf_confidence=[x for x in os.listdir(f'{ligand_fp}') if '_confidence' in x and x.endswith('.sdf') and x.startswith('rank')]
  num_predicts = len(sdf_confidence)

  # parse confidence score from sdf filename
  confidence_dict={}
  for rank in range(1,num_predicts+1):
    sdf_rank=[x for x in sdf_confidence if x.startswith(f'rank{rank}_confidence')][0]
    confidence_score=float(sdf_rank.split('confidence')[1].split('.sdf')[0])
    #print(f'reading confidence {confidence_score} from rank {rank} file: {sdf_rank}')
    confidence_dict[rank]=confidence_score

  # fetch color list
  rank_color = ligand_rank_color(cmap=cmap, data_list=confidence_dict.values())

  cmd.load(filename=protein_fp, object=protein_alias)

  for i in range(int(rank_num[0]), int(rank_num[1]) + 1):
    if os.path.exists(f'{ligand_fp}/rank{i}_reverseprocess.pdb'):
      # cmd.load( filename [,object [,state [,format [,finish [,discrete [,multiplex ]]]]]] )
      cmd.load(filename=f'{ligand_fp}/rank{i}_reverseprocess.pdb',
               object=f'{ligand_alias}_{i}',
               state=2,
               #discrete=1,
               multiplex=0,
               )

      # what a magic trick?!
      cmd.do(f'h_fix {ligand_alias}_{i}')
      print(f'color rank#{i} with color {rank_color[i-1]}')

      # set ligand backbone color
      cmd.set_color(f'color_ranked_{i}', rank_color[i-1])
      cmd.color(f'color_ranked_{i}', f'({ligand_alias}_{i})')
      util.cnc(f'({ligand_alias}_{i})', _self=cmd)

    else:
      # ignore if index of files is out of range
      print(f'file not found: {ligand_fp}/rank{i}_reverseprocess.pdb')


  cmd.set('surface_cavity_mode',4)
  cmd.set('transparency',0.5)
  cmd.set('surface_color','gray70')
  cmd.show('surface', protein_alias)
  cmd.hide('cartoon')

  cmd.set('stick_color','firebrick', protein_alias)


  cmd.center(f'{protein_alias} and not hetatm')

  os.makedirs('pse',exist_ok=True)
  cmd.save(f'pse/{protein_alias}_{ligand_alias}_{args.suffix}.pze')