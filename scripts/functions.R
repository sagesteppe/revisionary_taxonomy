
######################
###    reed_csv    ###
######################

reed_csv <- function(x){
  
  # This serves to read in and subset well behaved Darwin core data 
  # from SoRo note that the columns of interest can be readily changed. 
  # via selected cols of interests BY:
  # > WANTED_COLS <- paste(shQuote(names(TARGET_DF)), collapse=", ")
  
  wanted_cols <- c('id', 'institutionCode', 'ownerInstitutionCode', 
                   'catalogNumber', 'family', 'scientificName', 'taxonID',
                   'scientificNameAuthorship', 'identifiedBy', 'recordedBy', 
                   'associatedCollectors', 'recordNumber', 'eventDate', 'year',
                   'month', 'day', 'startDayOfYear', 'decimalLatitude', 
                   'decimalLongitude', 'geodeticDatum', 'georeferencedBy',
                   'georeferenceProtocol', 'georeferenceSources', 
                   'georeferenceVerificationStatus', 'georeferenceRemarks', 
                   'minimumElevationInMeters'
                   )
  
  # read in the data
  data <- read.csv(x, stringsAsFactors = F)
  
  # upon reading we do a few house keeping tasks, in theory
  # this will remove duplicates of the same collection in the
  # database but...
  output <- data %>% select(one_of(wanted_cols)) %>% 
    mutate(across(everything(), ~na_if(., ""))) %>% 
    drop_na(decimalLatitude, year)%>% 
    group_by(recordNumber, eventDate) %>% 
    distinct(.keep_all = T) %>% 
    distinct(recordNumber, .keep_all = T) %>% 
    filter(year >= 1800) %>% 
    mutate(recordedBy = str_remove(recordedBy, 'collectors: '))
  
  return(output)
}


###############################
###    population_finder    ###
###############################

# to this will need to be sure to not remove type populations
population_finder <- function(input, utm_epsg, combine_distance){
  
  # this functions seeks to identify herbarium specimens & 
  # other location based records which were taken from the
  # same group of individuals - for convenience a 'population'
  
  # input required is an sf object, in particuliar this was
  # developed on subsetted DarwinCore records. 
  # 'utm_epsg' refers to the code for a UTM zone - google can help
  # 'combine_distance' is how far this function will consider 
  # individual records to be paired based on them both being within
  # this distance. Pairs will then be subjected to graph building
  # to identify all individuals which share any linkages. 
  
  #################################################################
  
  # Each input is ensured to be in a numeric format.
  # the input crs is saved for reprojection of the output object.
  
  utm_epsg <- as.numeric(utm_epsg)
  combine_distance <- as.numeric(combine_distance)
  input <- st_as_sf(input)
  input_crs <- st_crs(input)

  # convert from a Geographic CRS to a Projected CRS for measurements with spdep
  object_utm <- st_transform(input, crs = utm_epsg) %>% 
    rowid_to_column()
  
  # identify all points within d2 of each other where 'dx' is in meters
  neighs <- spdep::dnearneigh(object_utm, d1 = 1, d2 = combine_distance)

  
  neighs <- as.list(neighs)
  names(neighs) <- c(1:length(neighs))
  neighs <- Map(cbind, neighs, rowNum = names(neighs))
  neighs <- map(neighs, as_tibble) 
  neighs <- bind_rows(neighs) # EXP
  #neighs <- data.table::rbindlist(neighs) # ? see if can become bind_rows(). 
  names(neighs) <- c('Neighbor', 'Row')
  neighs <- filter(neighs, Neighbor != 0) # maybe just bring back to base?

  # populations may proceed in an elongated fashion so that the 1st point
  # is only connected to the final points via intermediates. 
  # we will  create exhaustive networks of all joined groups of records
  
  my_edgelist <-igraph::graph.data.frame(neighs)
  groups <- lapply(igraph::groups(igraph::components(my_edgelist)), FUN =  data.frame)
  names(groups) <- c(1:length(groups))
  groups <- Map(cbind, groups, Group = names(groups))
  groups <- bind_rows(groups) # Was changed from rbindlist - no ill effect
  names(groups) <- c('Individual', 'Groups')
  groups$Individual <- as.numeric(as.character(groups$Individual))
  
  # now add the 'group' identifier back onto the projected sf/tibble
  object_utm <- right_join( groups,object_utm, by = c('Individual' = 'rowid')) 

  # We now select the most recent record to serve as the representative 
  # of each group identified by the process.  
  representatives <- object_utm %>% 
    drop_na(Groups) %>% 
    group_by(Groups) %>% 
    arrange(year) %>% 
    slice_max(year, with_ties = F)

  object_utm <-  object_utm %>% filter(is.na(Groups))
  output <- bind_rows(object_utm, representatives) %>% 
    st_as_sf()
  
  return(output)

}


#############################
###     closest_sites     ###
#############################


closest_sites <- function(input, target, congener, proximity_sites, random_sites){
  
  # the goal of this function is to find a certain number of the CLOSEST
  # sites to individuals of your target taxa, and to supplement these
  # with another number of sites which are randomly selected by a weight
  # to priotize closer sites, I usually run a 2/3 split between proximity
  # based and distance based. 
  
  # the idea being is the taxa are interfertile the more proximal non-target
  # populations should be closely related relative to the further ones. 
  
  # ==============================================
  
  
  target <- filter(input, scientificName == target)
  congener <- filter(input, scientificName == congener)
  
  # calculate the initial distances between ALL sites from the taxa
  distances <- st_distance(target, congener)
  rownames(distances) <- target$id
  colnames(distances) <- congener$id
  
  distances <- data.frame('distance' = sort(colSums(distances)/10)) %>% 
    rownames_to_column('id') %>% 
    mutate(id = as.numeric(id))
  
  # here we just find the closet sites to our target taxa
  closest <- distances[1:proximity_sites,]
  
  # now we randomly select some more - but weight closer one's still
  distances <- distances[(proximity_sites+1):nrow(distances),] %>% 
    slice_sample(n = random_sites, replace = F, weight_by = distance) %>% 
    rbind(., closest) %>% 
    mutate(Sample = 'Y')
  
  # Now combine the closest sites
  output <- left_join(input, distances, by = 'id') %>%
    st_as_sf() %>% 
    mutate(Sample = if_else(is.na(Sample), 'N', Sample)) %>% 
    dplyr::select(-distance) %>% 
    
  return(output)
}

#####################################
###     target_taxon_grp_fndr     ###
#####################################

target_taxon_grp_fndr <- function(input, exterior_sites, total_sites, mandatory_sites){
  
  # 'input' an sf/itbble dataset which contains the species of interest
  
  # 'id_col' = the column holding unambigious id's for these sites
  # 'exterior_sites' = the number of populations which form a geometric edge
  # convex hull around the species. 
  # 'total_sites' = the total number of sites you desire to sample from. 
  # 'mandatory_sites' = if applicable the ID of any sites you deem require be 
  # sampled, this should ALWAYS INCLUDE THE TYPE POPULATION 
  # 'unsamplable' = sites which you deem not possible to sample.

  distances <- function(x){
    distance_results <- sum(rowSums(as.matrix(st_distance(x)))) / length(x)
    distance_results <- distance_results/10
  }
  
  binomial <- input$Scientific_Name[1]
  # Identify all of the boundary sites. 
  species_boundary <- input %>% 
    st_union() %>% 
    st_convex_hull() %>% 
    st_cast("LINESTRING")
  boundary_groups <- unlist(st_intersects(species_boundary, input) )
  
  # combine the sites which will be kept
  pos_mand_sites <- which(input$id %in% mandatory_sites)
  all_kept_groups <- c(boundary_groups, pos_mand_sites)
  all_grps <- seq(1:nrow(input)) 
  grps_to_test <- setdiff(all_grps, 
                          all_kept_groups
                          )
  
  number_random_sites_needed <- total_sites - (as.numeric(length(boundary_groups))
  + as.numeric(length(mandatory_sites)))
  
  # determine the mean distance between all possible sites. 
  all_mean_dist <- rowSums(as.matrix(st_distance(input)))
  all_mean_dist <- sum(all_mean_dist)/length(all_mean_dist)
  
  # find all combinations of the needed remaining sites
  combinations <- combinat::combn(grps_to_test, 
                                  number_random_sites_needed) 
  all_kept_groups.df <- matrix(rep(all_kept_groups, ncol(combinations)),
                               nrow = length(all_kept_groups))
  combinations <- rbind(combinations, all_kept_groups.df)
  
  # repeatedly subset the sf dataframe to calculate the distances and
  # return a vector of ALL mean distances
  distances_vector <- rep("", times = ncol(combinations))
  for (i in 1:ncol(combinations)) {
    subset <- combinations[,i]
    distances_vector[i] <- distances(input[subset,]) # why no i in column?
  }
  names(distances_vector) <- 1:length(distances_vector)
  distances_vector <- sort(distances_vector, decreasing = T)[1:10]
  
  # now subset the combinations of sites to sample
  groups_combos <- combinations[,as.numeric(names(distances_vector))]
  
  # prettify results
  groups_combos <- as.data.frame(groups_combos)
  names(groups_combos) <- paste0(rep('Max_mean_dist', 10)
                                 , '_', seq(1, 10 ,1))
  
  # write out results
  binomial <- input$scientificName[1]
  name_write <- paste0(substr(binomial, 1,2), '.',
                       sapply(strsplit(binomial, ' '), `[`, 2))
  write.csv(groups_combos, 
            paste0(here(), '/data/processed/','Trgt_grp_fndr_', name_write ,'.csv'), 
            row.names = F)
  
  # return results
  return(groups_combos)

}

grps_most_dispersed <- function(x, max_oversample){
  
  # the purpose of this function is to identify which grps are most common 
  # in the maximum mean distances from each other. 
  
  
  # input is the output of 'target taxon grp finder'. Note row numbers
  # should not be present. 'max_oversample' refers to the maximum number of
  # over sample groups to return, this value not always achievable without 
  # altering parameters in the 'target taxon grp finder' function.
  
  # output of the function will be a vector in ascending 
  # order of the sites with the greatest mean distance. All unique groups from 
  # the analysis will be returned in case of need of oversampling
  
  number_grps <- length(unique(unlist(x)))
  max_grps <- ifelse(number_grps >= max_oversample, max_oversample, number_grps)
  
  groups <- as.numeric(names(sort(table(unlist(x)), decreasing =  T)[1:max_grps]))
  
  return(groups)
  
}




####################################
###         JOHN_TIGGER          ###
####################################

john_tigger <- function(x) {
  name <- x$NAME
  county_roads <- tigris::roads(state = x$STATEFP,
                                county = x$COUNTYFP, year = 2020)
  county_roads <- mutate(county_roads, NAME = name)
  
  return(county_roads)
}


###################################
####          PTS_CELL        #####
###################################

pts_cell <- function(x, utmzone, buf_dist, out_crs) {
  
  buffered <- x %>% 
    st_as_sf() %>% 
    st_transform(32613) %>% 
    st_buffer(14000) %>% 
    st_transform(4269) %>% 
    st_bbox() %>% 
    st_as_sfc() %>% 
    st_as_sf()
  
  return(buffered)
  
}

###################################
###        WINDOW_SEARCH       ####
#################################

window_search <- function(x){
  
  # The purpose of this function is to evaluate the buffered
  # 'windows'/panes/pixels which come from the 'pts_cell' function 
  # and try to minimize the number of windows required to fit around
  # all points.
  
  # as input it requires the output of an intersection between the output
  # of pts_cell and the original points which were fed into the pts_cell func.
  
  names(x) <- paste0('row', 1:nrow(x))
  
  # step 1 ######################################
  # IDENTIFY ALL SINGLETON WINDOWS 
  # then set aside for return later
  aLL <- lapply(x, length)
  singletons <- x[aLL <2] # identify ALL lone windows
  multiples  <- x[aLL>=2] # the inverse finds all multiples
  # these will be subject to testing for minimizing the number of windows
  ##############################################
  
  # step 2 #######################################
  # IDENTIFY ALL RECIPROCAL PAIRS - these will not under go network tests
  a <- map(x, data.frame) %>% 
    dplyr::bind_rows(.id = "window") %>% 
    rename('neighbor' = 2) %>% 
    group_by(window) %>% 
    mutate(row_num = 1:n()) %>% 
    pivot_wider(values_from = neighbor, names_from = row_num) %>% 
    unite("Pairs", 2:ncol(.), na.rm = TRUE, sep = '-') %>% 
    ungroup() %>% 
    group_by(Pairs) %>% 
    mutate(Grp_l = n()) %>% 
    mutate(identical_seq_l = str_count(Pairs, '-')+1) %>% 
    filter(Grp_l == identical_seq_l) %>% 
    filter(Grp_l > 1) 
  
  kept_pairs <- a %>% sample_n(1)
  pairs_remove <- a %>% filter(!window %in% kept_pairs$window)
  
  ############################################
  # part 3 subset original data set with removed points,
  kept <- paste(kept_pairs$window, collapse = "|")
  remove <- paste(pairs_remove$window, collapse = "|")
  
  groups <- multiples[str_detect(names(multiples), kept, negate = T)]
  groups <- groups[str_detect(names(groups), remove, negate = T)]
  
  # step 4 ########################################
  # identify groups with shared individuals
  b <- map(groups, data.frame) %>% 
    dplyr::bind_rows(.id = "window_id")
  colnames(b) <- c('window_id', 'target')
  b$target <- paste0('row',b$target)
  b <- b[b$window_id != b$target,]
  
  igraph_ob <- igraph::graph.data.frame(b, directed=FALSE)
  ident_grps <- lapply(igraph::groups(igraph::components(igraph_ob)),
                       FUN =  data.frame)
  
  names(ident_grps) <- c(1:length(ident_grps))
  ident_grps <- Map(cbind, ident_grps, Group = names(ident_grps))
  ident_grps <- bind_rows(ident_grps)
  colnames(ident_grps) <- c('rownum', 'Group')
  ident_grps <- tibble(ident_grps)
  
  vals <- lapply(groups, unlist)
  vals <- cbind(names(groups), vals)
  colnames(vals) <- c('rownum', 'vals')
  vals <- as_tibble(vals)
  vals$rownum <- as.character(unlist(vals$rownum))

  ident_grps <- left_join(tibble(ident_grps), tibble(vals), by = 'rownum')
  ident_grps <- drop_na(ident_grps)
  ident_grps <- tibble(ident_grps)	%>%
    rowwise() %>% 
    mutate(members_L = length(vals))	%>% 
    ungroup()
  
  # find the window with the most points in it 
  # and search for the most dissimilar sites below...
  colnames(ident_grps)[1] <- "rownum"
  network_grps <- split(ident_grps, ident_grps$Group)  
  
  group_partition <- function(x){
    
    first <- x %>% 
      slice_max(members_L, n =  1, with_ties = F) 
    
    # with this code we can find the values missing from the longest sequence
    vect <- sort(unique(c(unlist(x$vals))))
    func <- function(x){setdiff(vect, unlist(x))}
    missing_vals <- func(first$vals)
    
    # with this we test how many of the missing values 
    # are contained in each of the other windows. 
    my_intersect <- function(x){
      res <- intersect(unlist(x$vals),missing_vals) 
      res <- ifelse(suppressWarnings(res == missing_vals), 1, 0)
      res <- sum(res/length(res))
      return(res)
    }
    
    prop_rmn_cptr <- cbind(x, 
                           'NTRSCT' = apply(x, 1, my_intersect)
    )
    
    second <- prop_rmn_cptr %>% 
      slice_max(NTRSCT, n =  1, with_ties = F) %>% 
      bind_rows(first)
    
    grps_sub <- prop_rmn_cptr %>% 
      mutate('Keep' = rownum %in% second$rownum) %>% 
      filter(Keep == T) %>% 
      dplyr::select(rownum, Group)
    
    rm(vect, missing_vals, first, second)
    return(grps_sub)
  }
  
  maximal_grp_mmbrs <- bind_rows(lapply(network_grps, group_partition))
  
  grps_keep <- c(names(singletons),
                 maximal_grp_mmbrs$rownum,
                 kept_pairs$window
  )
  grps_keep <- paste(grps_keep, collapse = "|")
  output <- x[str_detect(names(x), grps_keep, negate = F)]
  return(output)
  
}
  


####################################
###        WINDOW CENTER         ###
####################################

window_center <- function(x, y){
  
  # this function is an immediate followup function for window search.
  # it deals with edge cases of pt dense areas which have redundant windows
  # left over from the above function. It required two inputs:
  
  # x = the original points we want to place in 'windows'/pane/pixels
  # y = the output from a second iteration run of pts_cell - the run which 
  # took the output of 'window_search'.
  
  pts <- x %>% mutate(PL_ID = 1:nrow(.))
  out <- st_intersection(pts, y)
  
  tp <- out %>% 
    ungroup() %>% 
    group_by(PL_ID) %>% 
    tally() %>% 
    filter(n >1)
  
  # 1 - Intersect all points to the remaining cells. 
  # 2 - determine which points are present in multiple windows
  # 3 - if a cell does NOT contain any novel points - drop it
  
  new_polys <-st_join(y, tp, left = F) %>% 
    distinct(POLY_ID, .keep_all = T) %>% 
    st_union() %>% 
    st_cast("POLYGON") %>% 
    st_as_sf() %>% 
    mutate(CAST_ID = 1:nrow(.))
  
  pt_grps <- st_intersection(x, new_polys) %>% 
    group_by(CAST_ID) 
  
  center_windows <- function(x){
    
    # this function does most the lifting of the centering.
    
    output <- x %>% 
      st_as_sf() %>% 
      st_union() %>% 
      st_convex_hull() %>% 
      st_transform(32613) %>% 
      st_buffer(1) %>% 
      st_cast('POLYGON') %>% 
      st_centroid() %>% 
      st_as_sf()
    
    return(output)
  }
  
  babadook <- split(pt_grps, ~ CAST_ID)  # MAKE THIS ACCEPT A NAMED COLUMN ... ?
  new_centers <- lapply(lapply(babadook, center_windows), pts_cell) %>% 
    bind_rows()
  
  test <- lapply(babadook, center_windows) %>% bind_rows()
  colnames(test) <- 'geometry'
  st_geometry(test) <- "geometry"
  sf <- test %>% 
    st_as_sf(sf_column_name = geometry) %>% 
    split(., 1:nrow(.))%>%
    purrr::map(pts_cell) %>% 
    bind_rows() %>% 
    st_as_sf()
  
  return(sf)
  
}

