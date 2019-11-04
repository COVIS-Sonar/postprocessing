#!/usr/bin/env python3
#
# I'm sure this is suuuper non-Pythoni
#
import glob
import argparse
import os.path

subst_dir = "substitute_functions/*.m"
subst_list = glob.glob( subst_dir )
subst_dict = {}
for s in subst_list:
    subst_dict[ os.path.basename(s) ] = s

cli = argparse.ArgumentParser()
cli.add_argument(
  "infiles",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*" )

# parse the command line
args = cli.parse_args()

def substitute(orig,subst_dict):
    bn = os.path.basename(orig)

    if bn in subst_dict:
        return "" #subst_dict[bn]
    else:
        return orig

outfiles = [substitute(a,subst_dict) for a in args.infiles]

print( ' '.join(outfiles) )
