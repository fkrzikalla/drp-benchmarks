drp-benchmarks
==============

Supplements to "Digital Rock Physics Benchmarks"

AUTHORS

Heiko Andrä
Jack Dvorkin
Matthias Kabel
Youngseuk Keehm
Minhui Lee
Junehee Han
Fabian Krzikalla (krz@stanford.edu)
Claudio Madonna
Mike Marsh
Nicolas Combaret
Tapan Mukerji
Sarah Ricker
Erik Saenger
Ratnanabha Sain
Nishank Saxena
Andreas Wiegmann
Erik Glatt
Xin Zhan

CONTENT

The natural samples are provided in raw grayscale tif format (Berea sandstone and Grosmont carbonate), and as segmented images. Segmented images are provided in raw unsigned char format (uint8), compressed using gzip.

./images/berea:
  segmented-kongju.mat     segmented-stanford.raw.gz  tif
  segmented-kongju.raw.gz  segmented-vsg.raw.gz
./images/fontainebleau:
  segmented-exxon.raw.gz
./images/grosmont:
  segmented-kongju.mat     segmented-su.raw.gz   tif
  segmented-kongju.raw.gz  segmented-vsg.raw.gz
  
For the generated shpere pack, we provide a digitized segmented image in compressed unsigned char format, as well as the list of sphere locations (final_loc_cn1.dat) and the size of the periodic box containing the 621 spherical objects. The matlab file discretizeSpherepack.m performs the discretization of the sphere pack geometry to a regular Cartesian grid. Six figures show the statistics and visualization of the pack.

./images/spherepack:
  Box_DIM.dat             figure03.eps  final_loc_cn1.dat
  discretizeSpherepack.m  figure04.eps  segmented-788x791x793.raw.gz
  figure01.eps            figure05.png
  figure02.eps            figure06.png

Laboratory results for the Mercury intrusion experiment are provided in open document spreadsheet (ods) format as well as in MS Excel format.

./labdata:
  bereasMercuryIntrusion.ods  bereasMercuryIntrusion.xls  figure02.pdf

The simulation results are to be found in simulations.ods / simulation.xls. We report the summary of the results as well as the digital rock physics simulation results performed on subvolumes.

./results:
  simulations.ods  simulations.xls
