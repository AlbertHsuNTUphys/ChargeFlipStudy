mkdir -p $2
python scripts/plotter.py -i   $1 -l 41500    -j test/analysis/ExYukawa/samples_2017.json  -o final_plotter.root -O $2 --normToData
cp test/index.php $2
