mkdir -p $2
python scripts/plotter.py -i $1 -l 41480.  -j test/analysis/ExYukawa/samples_2017_nanoaod.json  -o final_plotter.root -O $2 
cp test/index.php $2
