"""
eps2.tiff.py
Convert .eps files to .tiff files with the help of magick
Required: magick MUST BE INSTALLED and accessible from the command prompt
magick can be found at https://imagemagick.org/script/command-line-tools.php

Usage:
    eps2tiff.py [input eps dir] [output tiff dir]
"""


import sys
import os
import subprocess as sb

indir = sys.argv[1]
outdir = sys.argv[2]
files = os.listdir(indir)
files = [f for f in files if f.endswith('.eps') and not f.startswith('.')]
print(files)

for f in files:
    fn = f.split('.')[0]
    outfn = f"{fn}.tiff"
    cmd = ['magick', os.path.join(indir, f), os.path.join(outdir, outfn)]
    print(cmd)
    sb.run(cmd)

