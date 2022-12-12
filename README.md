# oncovirus

# how to use
  # place in same directory as files
  # use csv files
  # use an IDE/miniconda terminal
    # require downloads of numpy and pandas - conda install numpy
  # upon running files, will ask for input (rather than automatically running the analysis on a pre-specified file)
    # interaction!!
  # output files are pretty bare-bones, copy/paste into a "stats" excel worksheet and edit (add headers, sort)
    # will ask for names of output files - just give basic names, the point is to move the data in these files to one central file

# find motifs - identify_motif.py
  # bread and butter
  # input file should have data in the first column (i.e., not a column of 1, 2, 3, ...)
  # three options
    # normal: standard
    # peptide/protein full count: gives motifs per protein and per peptide as well as peptides per protein
    # abbreviated count: same as normal but each peptide is restricted to one instance of a particular motif
      # not that useful

# miscellaneous functions on csv sheets - identify_find.py
  # pretty much just an attempt to automate some of the more annoying functions
  # not that useful
  
# identify significantly enriched (not sure if also depleted) motifs between "populations" - fisher_motif.py
  # fisher's exact test
