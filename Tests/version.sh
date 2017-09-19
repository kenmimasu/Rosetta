cd '/Users/Ken/Work/Projects/Rosetta/dev/rosetta/Rosetta'
find . -name \*.pyc -delete
rm *.log
rm *.dat

cd '../'

tar zcf Rosetta-$1.tgz Rosetta

scp Rosetta-$1.tgz kenmimasu@login.hepforge.org:/hepforge/downloads/rosetta/