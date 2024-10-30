Hydrogen-bond Based Secondary Structure analysis program (HBSS).
This algorithm was designed to provide fine grained secondary structure (SS) classification based
on backbone hydrogen bond (Hb) patterns. HbSS is part of the SESCA package for CD spectrum
predictions or as a stand-alone python package for analysing protein structures.

License:  
HbSS and all its modules fall under the GNU general public license, version 3
(provided in the license.txt file in the main SESCA directory). 

This software is distributed in the hope that it will be useful but comes without any warranty.
The author does not accept responsibility to anyone for the consequences of using it, for whether
it serves any particular purpose, or that it works at all. No warranty is made about the software
or its performance. 




Contents:
The HbSS package contains four modules:
HBSS_main.py:   This is the main command line interface, allowing automated SS determination.
HBSS_prep.py:   Pre-processing module, analyses and re-label PDB files, calculates Hb patterns.
HBSS_basic.py:  Basic classification module, reads Hb pattern and classifies protein residues.
HBSS_extend.py: Advanced classification module, fine grained classification based Twist angles.

Note that each module is functional as a stand-alone python script, but can be imported into other
projects as a module. If no arguments are provided, each module prints its own usage and possible
command flags.




Installation:
HbSS is a python package which does not require compilation, and can run using both python 2.7
and python 3.4. It does not have any special dependencies besides the standard python modules
(sys, os, time and math).

To install HbSS:
1) Copy and unpack the zip file in your preferred directory (e.g.: '/home/usr/Programs/'), which
   creates a SESCA directory there
2) Enter the 'HBSS' subdirectory, and open the HBSS_main.py file with your favorite text editor
3) Modify the variable 'HBSS_dir' (in line 39) to point to the main HBSS directory
   (e.g.: HBSS_dir = "/home/usr/Programs/SESCA/HBSS")
4) Save the changes and close HBSS_main.py

It is also recommended to add the HBSS main directory to both your system path (e.g.  $PATH) and
your PYTHONPATH environmental variable for optimal use. If you wish to use HbSS as part of the
SESCA package for CD spectrum predictions, make sure that the SESCA_dir and HBSS_dir variables
are correctly configured in the HBSS_main.py module (found in SESCA/scripts).




Hydrogen bond definitions and default parameters (HBSS_prep):
HBSSS_prep.py analyses either individual protein structures in a standard PDB 1.0 format, or an
ensemble of several structures in a single PDB file, each inserted between lines starting with
'MODEL' and 'ENDMDL' respectively. Within each structure, HbSS searches for backbone hydrogen
bonds between amino acid residues.

By default HbSS assumes that model is full-atom model with hydrogens ( mode 1 in HbSS_prep).In
this case, a Hb between the residues i and j is identified if the distance between Hi and Oj is
smaller 2.3 Angstroms, and the angle Ni - Hi ---- Oj angle is larger than 135.0 degrees.
For this mode, backbone hydrogens are identified by the atom labels "H", "NH","H1","H2", and "H3",
backbone nitrogens should be labelled "N", and backbone oxygens as "O", "O1", "O2", "OT1", or "OT2"

For many PDB structures based on X-ray crystallography, hydrogen atoms are not represented
explicitly, and for these structures the presence of hydrogen bonds can be estimated from heavy
atom coordinates using operation mode 0 (or hb_type 0 in HBSS_main). In this case Hbs are estimated
between residues i and j if the Ni --- Oj distance is smaller than 3.3 Angstroms, and the angle
between CAi - Ni ---- Oj is larger than 100.0 degrees. CAi here refers to the C-alpha atom of
residue i, typically labelled "CA" in the PDB file. The labels for Ni and Oj are unchanged.

As a support function, HbSS prep can perform contact analysis between two atom groups (defined by)
the command flags @groupA and @groupB using operation mode 2. In this case all atoms in the defined
groups are considered potential contact partners if the distance to any atom form the other group is
within the cut-off distance of 3.5 Angstroms, and no angle limitations are set.

The described Hb definitions are defined between lines 215 - 245 in the HbSS_prep module, and can also
be customized using the command line flags @mode, @dist, and @angle.




HbSS  basic class definitions (HBSS_basic): 
Once PDB structure was analysed for hydrogen bond patterns, the amino acids of the protein are
classified into seven SS classes. In following we will denote a hydrogen bond between a donor
residue j and an acceptor residue i as "i -> j".

The basic classes for groups residues (henceforth segments) is defined as follows:
4-Helix segments (4H):
The 4-Helix class corresponds to the regular Alpha helix SS elements in other methods such as
DSSP, DISICL and STRIDE. In HbSS 4-Helix segments are identified based on local Hb patterns, a
segment of five residues (from i to i+4) is identified as a 4-Helix segment if there are Hb-s
between residues (i+4 -> i) and either (i+4 -> i+1) or (i+5 -> i+1)

3-Helix segments (3H):
The 3-Helix segments  correspond to the tighter 3/10 helical SS elements.
In HbSS, a segment of four residues (from i to i+3) is identified as a 3-Helix if there are Hb-s
between residues (i+3 -> i) and either (i+4 -> i+1) or (i+5 -> i+1). 

5-Helix sgements (5H):
The 5-Helix corresponds to  Pi-helical SS elements with an enlarged helix structure.
A segment of six residues (from i to i+5) are identified as a 5-Helix segment if there are Hb-s
between residues (i+5 -> i) and either (i+5 -> i+1) or (i+6 -> i+1).

Parallel beta strands (BSP):
Segments of hydrogen-bonded beta sheets are classified as parallel beta strands, if there are HB-s
between residues (j -> i) as well as (i+2 -> j), and |i-j| > 5. If these conditions are met, then
the segments (from i-1 to i+1) and (from j-1 to j+1) are classified as parallel beta strands.

Anti-parallel beta strands (BSA):
Segments of hydrogen-bonded beta sheets are classified as anti-parallel strands, if there are Hb-s
between residues (j -> i) as well as (i+2 -> j-2), and |i-j| > 6. The segments must also have Hb-s
between (i -> j) and (j-2 -> i+2). If these conditions are met, then the segments (from i to i+2)
and (from j-2 to j) are classified as anti-parallel beta strands.

Hydrogen bonded turns (TU):
Segments are classified as hydrogen bonded turns if they would not fit into any class from above,
but have a local Hb (i -> j) where  2 <= |i-j| >= 5. In this case, if i > j residues (from j to i)
are classified as a hydrogen bonded turn, otherwise the segment (from i to j) are classified as TU.

Unclassified (UNC):
All residues that are not classified into one of the six classes above remain unclassified.

Please note that if a residue would fit multiple classifications, it will assigned to the class
described first in this document (4H > 3H > 5H > BSP > BSA > TU > UNC).



Extended classification based on beta strand twist angles (HBSS_extend):
The HBSS_extend.py allows further extension of the beta strand classifications into left-handed,
relaxed, and right handed parallel or anti-parallel beta strands, increasing the number classes
to eleven.

Extending the classification involves the following steps:
1) The hydrogen bonded 3-residue segments are identified based on the basic HbSS.
   classification and the hydrogen bond output information.
2) For each segment pair the geometry of the backbone is extracted from the original PDB file.
3) The beta sheet twist angles are calculated for the 3-residue segment pairs.
4) The beta sheet residues are re-classified into one of the six extended beta strand classes. 

For parallel beta strands, the geometry of the segment is determined for the central residues
(i+1) and j. This requires the Bi vector between the midpoints of the bonds C(i)-N(i+1) and
C(i+1)-N(i+2), the Bj vector between the midpoints of C(j-1)-N(j) and C(j)-N(j+1), and the vector
Dji connecting CA(j) to CA(i+1). The twist angle is calculated according to Ho et al. as

Twist = signF * acos(Bi*Bj)/(|Bi|*|Bj|), where signF = sign[(Bi x Bj) * Dji]

For anti-parallel beta strands, the geometry of the segment is determined for the central residues
(i+1) and (j-1). This requires the Bi vector between the midpoints of the bonds C(i)-N(i+1) and
C(i+1)-N(i+2), the Bj vector between the midpoints of C(j-2)-N(j-1) and C(j-1)-N(j), and the vector
Dji connecting CA(j-1) to CA(i+1). The twist angle is calculated according to Ho et al. as

Twist = 180 - signF * acos(Bi*Bj)/(|Bi|*|Bj|), where signF = sign[(Bi x Bj) * Dji]

Once beta twist angles are determined, they are assigned to the central residues of both segments.
Since any beta strand may have up two neighbouring strands, the final twist angle of each residue is
calculated as the average of its assigned twist angles. If a beta strand residue has no assigned
twist angle of its own, twist angles from adjacent residues are averaged over instead.  

If its final twist angle is smaller than 3 degrees, residues are assigned to a left-handed beta
strand, if its twist angle is between 3 and 23 degrees, segments are assigned as a non-twisted
beta strand, and if the twist angle is larger than 23 degrees the segments are assigned as a
right-handed beta twist segment. The default boundaries for the beta twist angles can be set
separately for parallel and anti-parallel strands in HbSS_extend,  line 49



Output Information:
By default, the main output of an HbSS classification is a text file containing the time series of
each SS class. The time series contain the name of the class (in comments), the time code (model ID)
and residue number of each amino acid assigned to the class, terminated by an "&". Additionally, 
the basic output file contains the percentage of residues assigned to each class (in comments).
Further statistics on the number of classified residues and segments can be requested using the @det
flag. The name of the output file can specified using the @write flag, and detailed residue statis-
tics can also be written into a separate file using the @sum flag. Finally we note that, both
HBSS_prep and HbSS_extend can write output files containing Hb and twist angle summaries.


If you found HbSS useful during your work please cite the following article(s):

SESCA: Predicting Circular Dichroism Spectra from Protein Molecular Structures
Gabor Nagy, Maxim Igaev, Soren V. Hoffmann, Nykola C. Jones, Helmut Grubmuller
BioRxiv preprint, doi: https://doi.org/10.1101/279752


