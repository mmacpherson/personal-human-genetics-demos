stupid-human-genetics-tricks
============================

Learning about health, family, and ancestry from genotyping chip data.

To run the modules, you need some sample genotype data. Currently only
23andMe-format data is supported. The file may be uncompressed or
gzipped.

Sample public 23andMe-format genotypes may be obtained at
[https://opensnp.org/](OpenSNP).

Currently there are three programs available, that demonstrate the
basics of analyses that can be performed on genotype data.

## 1. Celiac Risk Report

*Example run, with genotype file named ``datafile``*

``> python celiac_risk.py datafile``

Looks up four markers reported to be associated with Celiac disease,
shows the odds ratio associated with each marker, and gives a combined
estimate of risk for Celiac disease.


## 2. Relatedness Test

*Example run, with two genotype files named ``datafile1`` and ``datafile2``*

``> python family_test.py datafile1 datafile2``

Prints an ASCII-art representation of genetic sharing along chromosome
1. Unrelated people will show no sharing; related people will show
   some or complete sharing depending on the nature of the
   relationship.


## 3. European Ancestry Classification

*Example run, with genotype file named ``datafile``*

``> python european_ancestry_inference.py datafile EUROSNP.genomeorder.frequencies.txt``

Evaluates the likelihoods of all six possible pairs of parents, where
the source populations are North West Europe, South East Europe, and
Ashkenazi Jewish. Meaning, it classifies the supplied genome according
to its ancestry, assuming that those are the only possible
populations.

Relies on the EUROSNP marker frequencies dataset from
[*Price et al* (2008)](http://genepath.med.harvard.edu/~reich/EUROSNP.htm).
