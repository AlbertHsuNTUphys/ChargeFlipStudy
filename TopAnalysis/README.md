# ChargeFlip Study

## Control Plot
The main code is src/ChargeFlip_nanoaod.cc. To run it, you need to first change the output directory in test/analysis/ExYukawa/Nano_steerAnalysis.sh.
Then to test the code,
```
sh test/analysis/ExYukawa/Nano_steerAnalysis.sh -o TEST
```
To run it on condor,
```
sh test/analysis/ExYukawa/Nano_steerAnalysis.sh -o SEL
```
To generate the control plot
```
sh plot_and_copy.sh <input_dir> <output_dir>
```


