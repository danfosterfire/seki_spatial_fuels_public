# README

## seki_cwd.csv

Tabular observations of coarse woody debris particles from each CWD transect as
described in the methods. 

Fields are:
  - plot_id: The unique plot identifier for the observation
  - transect_id: the unique transect identifier for the observation
  - az: The azimuth of the transect from plot center
  - location_m: The location in meters of the observation along the transect
  - diam_cm: the diameter in centimeters of the CWD particle intersecting the 
  transect
  - decay_class: the decay class of the CWD particle
  - comments: Any comments relating to the observation
  
## seki_fwddiams.csv

List of all diameters from FWD particles intersecting the 1m subtransects where 
FWD diameters were individually recorded.

Fields are:
  - plot_id: The unique plot identifier for the observation
  - transect_id: the unique transect identifier for the observation
  - az: The azimuth of the transect from plot center
  - diam_mm: The diameter in milimeters of the FWD particle
  
## seki_fwdtallies.csv

Tallies of FWD particles intersecting the 1m subtransects where FWD particles 
were tallied.

Fields are:
  - plot_id: The unique plot identifier for the observation
  - transect_id: the unique transect identifier for the observation
  - az: The azimuth of the transect from plot center
  - location_m: the location of the 1m subtransect along the main transect, recorded
  at the midpoint of the 1m subtransect
  - a1h: The tally of 1h fuel particles along the 1m subtransect which is parallel to 
  the main transect
  - a10h: The tally of 10h fuel particles along the 1m subtransect which is parallel to 
  the main transect
  - a100h: The tally of 100h fuel particles along the 1m subtransect which is parallel to 
  the main transect
  - b1h: The tally of 1h fuel particles along the 1m subtransect which is orthogonal to 
  the main transect
  - b10h: The tally of 10h fuel particles along the 1m subtransect which is orthogonal to 
  the main transect
  - b100h: The tally of 100h fuel particles along the 1m subtransect which is orthogonal to 
  the main transect
  - comments: Any comments relating to the observation
  
## seki_litterduff.csv

Observations of litter and duff depth at sample locations within each plot.

Fields are:
  - plot_id: The unique plot identifier for the observation
  - transect_id: the unique transect identifier for the observation
  - az: The azimuth of the transect from plot center
  - location_m: the location of the depth sample along the main transect
  - litter_cm: the depth of litter in centimeters
  - duff_cm: the depth of duff in centimeters
  - fuel_cm: the depth of the fuelbed in centimeters
  - comments: Any comments relating to the observation
  
## seki_metadata.csv

Metadata about each inventory plot.

Fields are:
  - plot_id: The unique plot identifier
  - transect_id: The unique transect identifier
  - az: The azimuth of the transect
  - slp_degrees: The slope along the transect from plot center to the point 30m 
  away at the end of the transect, measured in degrees
  - cwd: 0/1 whether coarse woody debris data were collected for the transect
  - fwd: 0/1 whether fine woody debris data were collected for the transect
  - litter_duff: 0/1 whether litter/duff depths were collected for the transect
  - veg: 0/1 whether understory vegetation data were collected for the transect
  - trees: 0/1 whether trees and saplings data were collected for the transect
  - comment: transect level comments
  - comment2: transect level comments
  - comment3: transect level comments

## seki_trees.csv

Treelist of trees and sapling falling within the belt transects.

Fields are:
  - plot_id: The unique plot identifier
  - tree_id: The unique tree identifier
  - ts: The "transect side" (two letters indicating cardinal directions, 
  first outward from plot center and then orthogonal to the first direction)
  - location_m: The distance from plot center along the first cardinal direction in 
  'ts'
  - distance_m: The distance from the main transect along the second cardinal direction 
  in 'ts'
  - status: the live/dead status of the tree
  - spp: The 4 letter species code of the tree
  - dbh_cm: The diameters in centimeters at breast height of the tree/sapling
  - height_m: The height in meters of the tree/sapling
  - htcb_m: The height to live crown base of the tree/sapling
  - decay_class: For snags, the decay class of the tree/sapling
  - comments: Any comments relating to the tree/sapling
  
## seki_veg.csv

Line-intercept observations of understory vegetation sampled along each 30m 
transect.

Fields are:
  - plot_id: the unique plot identifier
  - transect_id: the unique transect identifier
  - az: the azimuth of the transect from plot center
  - start_m: The distance along the transect (in meters) where the understory 
  vegetation initially intersects the transect
  - end_m: The distance along the transect (in meters) where the understory 
  vegetation stops intersecting the transect
  - spp: The species of the understory vegetation as a four letter code
  - status: the live/dead status of the understory vegetation
  - height_m: The average height (in meters) of the understory vegetation over 
  its intersection with the transect
  - comments: any comments relating to the observation
  
  
