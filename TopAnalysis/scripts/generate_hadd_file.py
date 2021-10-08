import os
import sys
import optparse
import subprocess
import json

if __name__=='__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-f','--json', dest='json', help='json file', default=None, type='string')
    (args,opt) = parser.parse_args()


    with open(args.json, "r") as f:
        jsons = json.load(f)
        f.close()
    with open("./nano_hadd_list.sh", "w") as shell:
        shell.write("mkdir $2\n")
        for dataset in jsons:
            name = dataset['dir']
            print(name)
            shell.write("hadd $2/%s.root $1/%s_*.root\n"%(name,name))


