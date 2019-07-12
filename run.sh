source activate protac

 'This script runs the process of generating predictions for a PROTAC mediated
  ternary structure. The protein structure files, pdb_1 and pdb_ are docked
  using ZDOCK. Each of these files contains their respective ligands which are
  named, pdb_1_lig and pdb_2_lig. ZDOCK does not use hydrogens during docking,
  but prior to generation of output strucutres, filenames are provided which
  contain pdb_1 and pdb_2 with hydrogens added. This allows the complexes to
  be used with other software or progams that generally expect hydrogen atoms.
  The generated ZDOCK output structures are filtered using filter_zdock.py
  based on the provided  distance, dist, and by having >=100 angstroms^2 of
  hydrophobic solvent-accessible surface area buried. For structures which pass
  this filter, an attempt is made to generate a PROTAC conformer, based on the
  provided protac_smi arguement, using gen_conf.py which connects the binding
  sites of the docked protein-protein complex. These conformers are generated
  using a restrained  embedding that uses the maximum common substructures of 
  the ligands present  in each binding site, pdb_1_lig and pdb_2_lig, as atomic
  restraints. The complexes for which a PROTAC conformer was successfully
  generated can then compared a reference pdb structure if one was provided.


  Prior to running this script the following files are required:
	A zdock/ directory containing the zdock binary and associated files
	
  	The two .pdb files corresponding to each protein in the ternary complex
     with and without hydrogens (4 .pdb files total) in the zdock/ directory
	 -These files should ideally be named $name1_receptor.pdb,
	  $name1_receptor_H.pdb, $name2_ligand.pdb, and $name2_ligand_H.pdb 
	 -"receptor" and "ligand" here refer to how they are treated by ZDOCK 
        during the protein-protein docking process
	 -Using a different file naming scheme may require editing this script
	
	The names of the ligand present in each protein structure binding site
     need to be known and provided (2 strings)
     
	The smiles string of the PROTAC molecule of interest (1 string)
	
	An approximate distance of the length of the fully elongated PROTAC
	molecule (1 float/int)
	
	The zrank binary in the top level directory	
	
	Optional: The name of a reference .pdb file to compare the predicted
 			structures to.


  Example command:
	./run.sh BRD4_receptor.pdb CRBN_ligand.pdb JQ1 THL \
	"c1cc(Cl)ccc1C(c(c(C)c(s2)C)c2n(c34)c(nn4)C)=N[C@H]3CC(=O)NCCCCCCCCNC(=O)COc(ccc5)c(c56)C(=O)N(C6=O)[C@@H]7C(=O)NC(=O)CC7" \
	35 ./prep/6BN7_loopfill_H.pdb


  Note: This file contains many bash commands which are only necessary to 
	   connect outputs and inputs from the different python scripts. ZDOCK and
        the .py files are doing the actual important calculations. However,
	   the gen_conf.py script is run in parallel using the xargs bash command.
	   Running gen_conf.py in parallel is not strictly necessary, but can
	   greatly reduce runtime if there are a large number of complexes to 
        attempt to generate conformers for.
'

#Variables to store the positional command line arguements
pdb_1=`echo $1 | cut -f1 -d '.' | cut -f1 -d '_'`
pdb_2=`echo $2 | cut -f1 -d '.' | cut -f1 -d '_'`
pdb_1_lig=$3
pdb_2_lig=$4
protac_smi=$5
dist=$6
ref_pdb=$7

#Perform protein-protein docking
cd zdock
./mark_sur $1 "$pdb_1"_m.pdb
./mark_sur $2 "$pdb_2"_m.pdb
./zdock -R "$pdb_1"_m.pdb -L "$pdb_2"_m.pdb -o zdock.out -N 2000

#Replace protein structure file name for those with hydrogens in output file
# before using docking results to generate complexes
sed -i -e "s/"$pdb_1"_m.pdb/$pdb_1"_receptor_H".pdb/g" zdock.out
sed -i -e "s/"$pdb_2"_m.pdb/$pdb_2"_ligand_H".pdb/g" zdock.out

./create.pl zdock.out
mkdir out
mv complex.*.pdb out/
mkdir out_filter

#Filter complexes based on center of mass distances between ligands and
# the burial of hydrophobic solvent-accessible surface area
for f in ./out/*; do 
	pdb=`echo $f | xargs -n 1 basename`;
	python ../py/filter_zdock.py -p $pdb -a $pdb_1"_receptor_H".pdb \
	-b $pdb_2"_ligand_H".pdb -l $pdb_1_lig $pdb_2_lig -d $dist -i out/ \
	-o out_filter/ >> filter_zdock.log;
done

sed -i '/ PyMOL not running, entering library mode (experimental)/d' \
./filter_zdock.log


#For each structure, attempt to place a PROTAC molecule conformer to form the
# ternary complex 
cd ..
mkdir confs
cd confs/
for f in ../zdock/out_filter/*;do
	d=`echo $f | cut -f4 -d '.'`;
	mkdir $d;
	grep $pdb_1_lig $f > $d/$pdb_1_lig.pdb;
	grep $pdb_2_lig $f > $d/$pdb_2_lig.pdb;
done
cd ..

for d in confs/*; do 
	echo "source activate protac; python ./py/gen_conf.py \
	-p \'$protac_smi\' -l $pdb_1_lig.pdb -t $pdb_2_lig.pdb \
	-n 1 -c 0 -i $d -o $d >> gen_conf.log" >> command.txt;
done

cores=$(fgrep -c processor /proc/cpuinfo)
xargs --arg-file=command.txt \
      --max-procs=$cores  \
      --replace \
      --verbose \
      /bin/sh -c "{}"

cat gen_conf.log | grep complex > gen_conf.out


#Score complexes using ZRANK
cat gen_conf.out | while read line; do
	echo ./zdock/out_filter/$line >> gen_conf;
done 
./zrank gen_conf
rm  gen_conf


#Calculate RMSD to a reference pdb if provided
if [ -n "${ref_pdb+set}" ]; then
	python ./py/compare_ternary.py -p gen_conf.out \
	-i zdock/out_filter/ -r $ref_pdb
fi
