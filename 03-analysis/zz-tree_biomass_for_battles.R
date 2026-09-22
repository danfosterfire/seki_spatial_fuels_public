
table9 = read.csv(here::here('02-data',
                             '00-source',
                             'gill',
                             'table9.csv'))

lockhart_table_3 = 
  read.csv(here::here('02-data',
                      '00-source',
                      'lockhart',
                      'lockhart_table_3.csv'))



ref_species = 
  read.csv(here::here('02-data',
                      '00-source',
                      'fia',
                      'REF_SPECIES.csv')) %>%
  as_tibble() %>%
  
  select(GENUS, SPECIES,
         SPECIES_SYMBOL, 
         JENKINS_TOTAL_B1, JENKINS_TOTAL_B2,
         JENKINS_STEM_WOOD_RATIO_B1, JENKINS_STEM_WOOD_RATIO_B2,
         JENKINS_STEM_BARK_RATIO_B1, JENKINS_STEM_BARK_RATIO_B2,
         JENKINS_FOLIAGE_RATIO_B1, JENKINS_FOLIAGE_RATIO_B2,
         JENKINS_SAPLING_ADJUSTMENT, 
         MC_PCT_GREEN_WOOD,
         MC_PCT_GREEN_BARK,
         STANDING_DEAD_DECAY_RATIO1,
         STANDING_DEAD_DECAY_RATIO2,
         STANDING_DEAD_DECAY_RATIO3,
         STANDING_DEAD_DECAY_RATIO4,
         STANDING_DEAD_DECAY_RATIO5)

tree_sapling_biomass = 
  readRDS(here::here('PATHTOYOURTREELIST')) %>%
  # estimate the crown radius of each tree using the results of gill et al. 2000
  # first join in the regression coefficients
  left_join(table9, by = c('spp' = 'Species')) %>%
  
  # fill in unmatched values with the average of the hardwood species from the 
  # lockhart paper; which are different hardwood species but probably better than 
  # a generic conifer coefficient
  mutate(b0 = 
           ifelse(is.na(b0), 
                  lockhart_table_3[lockhart_table_3$species=='average','a'],
                  b0),
         b1 = 
           ifelse(is.na(b1), 
                  lockhart_table_3[lockhart_table_3$species=='average','b'], 
                  b1)) %>%
  
  # now, use the DBH and regression coefficients to estimate the crown radius
  mutate(crownrad_m = b0+(b1*dbh_cm)) %>%
  select(-b0, -b1) %>%
  
  # map the species codes we used to the official species codes used by the FIA
  # table; they don't appear to have generic coefficients to use for unknown 
  # species, so I'm just going to use the ABCO coefficients when species 
  # is unknown, based on the overstory composition data from seki
  left_join(data.frame(spp = c('ABCO', 'ACMA', 'ALRH', 'CADE',
                               'CONU', 'PILA', 'PIPO', 'QUCH', 'QUKE', 'TOCA',
                               'UNK'),
                       spp4 = c('ABCO', 'ACMA3', 'ALRH2', 'CADE27', 'CONU4',
                                'PILA', 'PIPO', 'QUCH2', 'QUKE', 'TOCA',
                                'ABCO')) %>%
              mutate(spp = as.character(spp),
                     spp4 = as.character(spp4))) %>%
  
  # join in the biomass coefficients
  left_join(ref_species %>%
              select(-GENUS, -SPECIES),
            by = c('spp4' = 'SPECIES_SYMBOL')) %>%
  
  # estimate the live biomass using equations from
  # PNW. (2015). A Data Dictionary and User Guide for the PNW-FIADB database. 
  # Pacific Northwest Research Station. 
  # https://www.fs.fed.us/pnw/rma/fia-topics/documentation/documents/PNW_FIADB_P2_Manual_2014.pdf
  # appendix J
  mutate(
    DIA = dbh_cm / 2.54,
    agbiomass_lb = 
      exp(JENKINS_TOTAL_B1+JENKINS_TOTAL_B2*log(DIA*2.54))*2.206,
    stem_ratio = 
      exp(JENKINS_STEM_WOOD_RATIO_B1+JENKINS_STEM_WOOD_RATIO_B2/(DIA*2.54)),
    bark_ratio = 
      exp(JENKINS_STEM_BARK_RATIO_B1+JENKINS_STEM_BARK_RATIO_B2/(DIA*2.54)),
    foliage_ratio = 
      exp(JENKINS_FOLIAGE_RATIO_B1+JENKINS_FOLIAGE_RATIO_B2/(DIA*2.54))) %>%
  
  mutate(
    # the foliage ratio equation breaks down for very small trees, giving 
    # ratios well over 1; cap it such that bark+stem+foliage=1
    foliage_ratio = 
      ifelse(bark_ratio+foliage_ratio+stem_ratio>1,
             1-bark_ratio-stem_ratio,
             foliage_ratio),
    stembiomass_lb = 
      agbiomass_lb*stem_ratio,
    barkbiomass_lb = 
      agbiomass_lb*bark_ratio,
    foliagebiomass_lb = 
      agbiomass_lb*foliage_ratio,
    branchbiomass_lb = 
      agbiomass_lb-stembiomass_lb-barkbiomass_lb-foliagebiomass_lb,
    
  ) %>%

  # adjust for snags
  mutate(
    
    # stem is woody, gets woody decay
    stembiomass_lb = 
      ifelse(status=='D',
             ifelse(decay_class=='1',
                    stembiomass_lb*STANDING_DEAD_DECAY_RATIO1,
                    ifelse(decay_class=='2',
                           stembiomass_lb*STANDING_DEAD_DECAY_RATIO2,
                           ifelse(decay_class=='3',
                                  stembiomass_lb*STANDING_DEAD_DECAY_RATIO3,
                                  ifelse(decay_class=='4',
                                         stembiomass_lb*STANDING_DEAD_DECAY_RATIO4,
                                         ifelse(decay_class=='5',
                                                stembiomass_lb*STANDING_DEAD_DECAY_RATIO5,
                                                stembiomass_lb))))),
             stembiomass_lb),
    
    # assume bark decays at the same rate as the wood (not great but also
    # not that important, and snag decay class provides little info about
    # how much bark is remaining)
    barkbiomass_lb = 
      ifelse(status=='D',
             ifelse(decay_class=='1',
                    barkbiomass_lb*STANDING_DEAD_DECAY_RATIO1,
                    ifelse(decay_class=='2',
                           barkbiomass_lb*STANDING_DEAD_DECAY_RATIO2,
                           ifelse(decay_class=='3',
                                  barkbiomass_lb*STANDING_DEAD_DECAY_RATIO3,
                                  ifelse(decay_class=='4',
                                         barkbiomass_lb*STANDING_DEAD_DECAY_RATIO4,
                                         ifelse(decay_class=='5',
                                                barkbiomass_lb*STANDING_DEAD_DECAY_RATIO5,
                                                barkbiomass_lb))))),
             barkbiomass_lb),
    
    # branches are woody, get woody decay
    branchbiomass_lb = 
      ifelse(status=='D',
             ifelse(decay_class=='1',
                    branchbiomass_lb*STANDING_DEAD_DECAY_RATIO1,
                    ifelse(decay_class=='2',
                           branchbiomass_lb*STANDING_DEAD_DECAY_RATIO2,
                           ifelse(decay_class=='3',
                                  branchbiomass_lb*STANDING_DEAD_DECAY_RATIO3,
                                  ifelse(decay_class=='4',
                                         branchbiomass_lb*STANDING_DEAD_DECAY_RATIO4,
                                         ifelse(decay_class=='5',
                                                branchbiomass_lb*STANDING_DEAD_DECAY_RATIO5,
                                                branchbiomass_lb))))),
             branchbiomass_lb),
    
    # assume all snags have no foliage
    foliagebiomass_lb = 
      ifelse(status=='D',
             0,
             foliagebiomass_lb)
  ) %>%
  
  # recalculate biomass aggregatsions after doing the decay adjustments
  mutate(
    agbiomass_lb = 
      stembiomass_lb + barkbiomass_lb + branchbiomass_lb + foliagebiomass_lb,
    
    # lump the bark in with the stem
    bolebiomass_lb = stembiomass_lb+barkbiomass_lb,
    
    crownbiomass_lb = branchbiomass_lb+foliagebiomass_lb) %>%
  
  # convert to metric for the values we intend to keep
  mutate(
    agbiomass_kg = agbiomass_lb *0.453592,
    bolebiomass_kg = bolebiomass_lb* 0.453592,
    crownbiomass_kg = crownbiomass_lb * 0.453592,
    foliagebiomass_kg = foliagebiomass_lb * 0.453592
  ) %>%
  
  # ditch all the intermediate columns
  select(
    group_id, plot_id,  x_coord, y_coord,
    status, spp, dbh_cm, height_m, htcb_m, decay_class,
    crownrad_m, agbiomass_kg, bolebiomass_kg, foliagebiomass_kg, crownbiomass_kg) %>%

  # calculate some standard summary statistics values
  mutate(
    tph = ifelse(dbh_cm>=11.4,
                 1/0.05,
                 1/((2*30+(2*28))/10000)),
    ba_m2 = pi*((dbh_cm/2/100)**2),
    ba_m2ha = ba_m2 * tph,
    agbiomass_mgha = 
      (agbiomass_kg/1000)*tph
  )