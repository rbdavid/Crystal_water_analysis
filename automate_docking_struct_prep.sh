
# step 1: run dx analysis on available crystal structures
echo "python3 dx_analysis.py dx_analysis.config IO.py"
python3 dx_analysis.py dx_analysis.config IO.py

AMBERHOME="/home/rbdavid/Apps/amber20"

FILES=structures/6WX4_chain_D_aligned_6WX4.pdb
selection='not type H ZN and (not segid B or resname HOH)'
mobile_landmark='name CA and resid 63:72 78:91 110:121 130:141 144:155 165:175'
alignment_landmark='name CA and resid 63:72 78:91 110:121 130:141 144:155 165:175'
for f in $FILES
do
	temp=${f##*/}
	base_file_name=${temp%%.*}	# rmvs .pdb extension characters
	echo $base_file_name
	mkdir $base_file_name

	# step 2: prepare the pre-aligned crystal structure for quick AMBER minimization of hydrogens
	sed -e "s?AAA?$f?g" -e "s?BBB?$selection?g" -e "s?CCC?$base_file_name/$base_file_name.pdb?g" prep_crystal.config > $base_file_name/prep_crystal.config
	echo "python3 prep_crystal.py prep_crystal.config IO.py"
	python3 prep_crystal.py $base_file_name/prep_crystal.config IO.py
	cd $base_file_name
	sed -e '/^CONECT/d' $base_file_name.pdb > temp.pdb	# -e '/^.*ZN.*/i TER' -e '/^.*ZN.*/a TER' 
	mv temp.pdb $base_file_name.pdb

	# step 3: send the pdb through tleap to prepare prmtop and inpcrd files
	sed -e "s/AAA/$base_file_name.pdb/g" ../tleap.in > temp.in
	echo "tleap -s -f temp.in > tleap.out"
	tleap -s -f temp.in > tleap.out

	# step 4: minimize the positions of all hydrogens (protein and waters) using a quick sander minimization
	echo "time $AMBERHOME/bin/sander -O -i ../min.in -o min.out -r min.rst -x min.nc -c protonated.inpcrd -p protonated.prmtop -ref protonated.inpcrd -inf min.mdinf"
	time $AMBERHOME/bin/sander -O -i ../min.in -o min.out -r min.rst -x min.nc -c protonated.inpcrd -p protonated.prmtop -ref protonated.inpcrd -inf min.mdinf

	# step 5: take last frame of minimization, analyze the positions of crystal waters (oxygen atoms), and see if those atoms pass the boolean tests defined in dx_booleans list of lists
	sed -e "s/AAA/$base_file_name/g" -e "s?BBB?../$f?g" -e "s/CCC/$mobile_landmark/g" -e "s/DDD/$alignment_landmark/g" ../prep_structures.config > prep_structures.config
	cp ../IO.py .
	echo "python3 prep_structures.py prep_structures.config IO.py"
	python3 ../prep_structures.py prep_structures.config IO.py
	cd ../
done

