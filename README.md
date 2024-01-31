# oncovirus

## how to use
- place in same directory as files
- use csv files
- use an IDE/miniconda terminal
- require downloads of numpy and pandas - conda install numpy
- upon running files, will ask for input (rather than automatically running the analysis on a pre-specified file)
- output files are pretty bare-bones, copy/paste into a "stats" excel worksheet and edit (add headers, sort)
- will ask for names of output files - just give basic names, the point is to move the data in these files to one central file

## identify_motif.py - find motifs
- bread and butter
- input file should have data in the first column (i.e., not a column of 1, 2, 3, ...)
- three options
  - normal: standard
  - peptide/protein full count: gives motifs per protein and per peptide as well as peptides per protein
  - abbreviated count: same as normal but each peptide is restricted to one instance of a particular motif
    - not that useful
    
## identify_find.py - miscellaneous functions on csv sheets
- pretty much just an attempt to automate some of the more annoying functions
- not that useful

## fisher_motif.py - identify significantly enriched (not sure if also depleted) motifs between "populations"
- fisher's exact test
