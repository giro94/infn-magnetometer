Usage:
./bin/fullMC <radius [0, 1]> <space model> <kicker transient hist name> <beam hist name>

Radius 0 is at the magic radius, 1 is at 17.5 mm

Space model can be:
  flat: flat space distribution, no kick weighting
  x2: x^2 distribution, y flat
  x4: x^4 distribution, y flat
  x4y2: x^4, y^2 distribution 
  x2m2: x^2 distrib with -2mm shift (for uncertainty)
  x2p2: x^2 distrib with +2mm shift (for uncertainty)
  x4m2: x^4 distrib with -2mm shift (for uncertainty)
  x4p2: x^4 distrib with +2mm shift (for uncertainty)
  hist: uses 2d histogram distribution in the file KickerSpaceModel_1.root named h2Space
  
Kicker transient hist name follows the naming Paolo gave (see INFN_Umass_hd.root file):
  hist_name           Vibr    Smooth
  h1_kick1_R0_ra:     yes,    yes
  h1_kick1_R0_p:      no,     no
  h1_kick1_R0_p_ra:   no,     yes
  h1_kick1_R0:        yes,    no
  h1_kick1_R1_ra:     no,     yes
  h1_kick1_R1:        no,     no
  
Beam hist name is the name of the beam distribution hist in the file BeamDistrib/beam_dists.root:
  noRF 
  xRF
  xyRF
  
The output file is written following the parameters given:
full_r<radius>_m<spacemodel>_K-<kickertransient>_B-<beamdist>.root 
