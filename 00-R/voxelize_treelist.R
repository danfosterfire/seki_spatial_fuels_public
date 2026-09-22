voxelize_trees_and_saplings = 
  
  function(treelist, voxel_dimensions){
    
    # for testing
     treelist = tree_sapling_biomass
     #must have plot_id, dbh_m, height_m, spp,
    # htcb_m, bolebiomass_kg, crownbiomass_kg, x_coord, y_coord, crownrad_m
    voxel_dimensions = 
      list(x_coords = seq(from = 0, to = 30, by = 0.5),
           y_coords = seq(from = 0, to = 30, by = 0.5),
           z_coords = seq(from = 0, to = 75, by = 0.5))
    #### bole radii parameters #################################################
    # for each tree, get the parameters a and b in the equation 
    # radius = a + b*height, 
    # for the bole: using the points (dbh_m, 1.37) and (0, height_m)
    # dbh_m = a + b*1.37               0 = a + b*height_m
    # a = dbh_m - b*1.37               0 = dbh_m - b*1.37 + b*height_m
    #                                  0 = dbh_m + b *(height_m - 1.37)
    #                                  b = -dbh_m / (height_m - 1.37)
    # a = dbh_m - (-dbh_m / (height_m-1.37))*1.37
    # 
    # for the crown: using the points (crownrad_m, htcb_m) and (0, height_m)
    # crownrad_m = a + b*htcb_m         0 = a+b*height_m
    # a = crownrad_m - b*htcb_m
    #                                   0 = crownrad_m - b*htcb_m + b*height_m
    #                                   0 = crownrad_m + b(height_m - htcb_m)
    #                                   -crownrad_m / (height_m - htcb)m) = b
    # a = crownard_m - (-crownrad_m/(height_m-htcb_m))*htcb_m
    
    treelist = 
      treelist %>%
      mutate(radius_bh_m = dbh_cm/200,
             a_bole = radius_bh_m - (-radius_bh_m / (height_m-1.37))*1.37,
             b_bole = -(radius_bh_m / (height_m-1.37)),
             a_crown = crownrad_m + htcb_m*(crownrad_m/(height_m-htcb_m)),
             b_crown = -(crownrad_m/(height_m-htcb_m)))
    
  
    # get volume of bole and crown, get density of bole and crown
    treelist = 
      treelist %>%
      mutate(bolevol_m3 = pi*((a_bole+(b_bole*0))**2) * (height_m/3),
             crownvol_m3 = 
               # cone for conifers, cylinder for hardwoods
               ifelse(is.element(spp, c('ABCO', 'CADE', 'PILA', 'PIPO', 'PSME')),
                      # cone for conifers
                      pi*(crownrad_m**2)*((height_m-htcb_m)/3),
                      # cylinder for hardwoods
                      pi*(crownrad_m**2)*((height_m-htcb_m))),
             boledensity_kgm3 = bolebiomass_kg / bolevol_m3,
             crowndensity_kgm3 = crownbiomass_kg / crownvol_m3) %>%
        # snags have NA crowndensity, which causes problems below. we want to treat 
      # snags as having 0 crown density
      mutate(crowndensity_kgm3 = 
               ifelse(is.na(crowndensity_kgm3),
                      0,
                      crowndensity_kgm3))
    
    ordered_plots = 
      unique(treelist$plot_id)[order(as.numeric(unique(treelist$plot_id)))]

    # for testing
    #ordered_plots = 1
    # more memory efficient version?
    bulkdensity_array = 
  
      array(dim = c(length(voxel_dimensions$x_coords),
                    length(voxel_dimensions$y_coords),
                    length(voxel_dimensions$z_coords),
                    length(ordered_plots)),
            
            dimnames = 
              list('x' = voxel_dimensions$x_coords,
                   'y' = voxel_dimensions$y_coords,
                   'z' = voxel_dimensions$z_coords,
                   'plot_id' = ordered_plots),
            
            data = 
              
              # for every plot
              sapply(X = 1:length(ordered_plots),
                     FUN = function(i){
                       
                       # get the trees on this plot
                       trees = which(treelist$plot_id==ordered_plots[i])
                       
                       spp = treelist$spp[trees]
                       
                       # get the vector of coordinates for each tree
                       x_tree = treelist$x_coord[trees]
                       y_tree = treelist$y_coord[trees]
                       
                       # get the htcb for each tree
                       htcb = treelist$htcb_m[trees]
                       height_m = treelist$height_m[trees]
                       crownrad_m = treelist$crownrad_m[trees]
                       
                       # get the bole radius coefficients for each tree
                       a_bole = treelist$a_bole[trees]
                       b_bole = treelist$b_bole[trees]
                       boledensity_kgm3 = treelist$boledensity_kgm3[trees]
                       a_crown = treelist$a_crown[trees]
                       b_crown = treelist$b_crown[trees]
                       crowndensity_kgm3 = treelist$crowndensity_kgm3[trees]
                       
                       # for every z
                       sapply(X = 1:length(voxel_dimensions$z_coords),
                              FUN = function(j){
                                
                                # get the z coordinate
                                z = voxel_dimensions$z_coords[j]
                                
                                # calculate the bole radius of each tree at this z,
                                # eqn goes negative above height of each tree, 
                                # so multiply by 0 if z is greater than each trees height
                                bole_radii = 
                                  (a_bole+(b_bole*z))*(z <= height_m)
                                  
                                # calculate the crown radius of each tree at this z
                                # multiply by 0 if z is greater than the tree's height,
                                # or below the trees htcb
                                crown_radii = 
                                  ifelse(is.element(spp, 
                                                    c('ABCO', 'CADE', 'PILA', 'PIPO')),
                                         # cone for conifers
                                         (a_crown+(b_crown*z))*(z>=htcb & z<=height_m),
                                         # cylinder for hardwoods
                                         crownrad_m * (z>=htcb & z<=height_m))
                                
                                # set NAs to 0 (for snags)
                                crown_radii[is.na(crown_radii)] = 0
                                
                                # for every y
                                sapply(X = 1:length(voxel_dimensions$y_coords),
                                      FUN = function(k){
                                        
                                        # get the y coordinate
                                        y = voxel_dimensions$y_coords[k]
                                        
                                        # for every x
                                        sapply(X = 1:length(voxel_dimensions$x_coords),
                                               FUN = function(l){
                                                 
                                                 # get the x coordinate
                                                 x = voxel_dimensions$x_coords[l]
                                                 
                                                 # get the distance between this 
                                                 # point and each tree
                                                 horiz_distance = 
                                                   sqrt((x-x_tree)**2+
                                                          (y-y_tree)**2)
                                                 
                                                 # get the (within plot) indices
                                                 # of the trees which this point 
                                                 # is within the crown
                                                 within_crown = 
                                                   (horiz_distance < crown_radii)
                                                 
                                                 within_bole = 
                                                   (horiz_distance < bole_radii)
                                                 
                                                 # get the contribution of each 
                                                 # tree to the total bulk density
                                                 bulkdensity_kgm3 = 
                                                   (within_bole*boledensity_kgm3)+
                                                   (within_crown*crowndensity_kgm3)
                                                 
                                                 # sum all the tree's contributions
                                                 # to get the overall bulk density 
                                                 # for this voxel
                                                 return(sum(bulkdensity_kgm3))
                                                 
                                               })
                                      })
                                
                                
                                
                              })
                       
                       
                     }))
        
    return(bulkdensity_array)
  }