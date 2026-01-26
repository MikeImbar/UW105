library(gtools)
library(matlib)
library(geomorph)
library(Morpho)
library(princurve)
library(tripack)
library(stringr)
library(geometry)
library(randomcoloR)
library(Rvcg)

# UM1/UM2/UM3/LM1/LM2/LM3
position <- c("UP", "UP3", "UP4")

# Choose an analysis type. 
# Options: "EDJ CEJ", "EDJ", "CEJ"

analysis <- "EDJ CEJ"

# Remove the landmark type and file extension from the landmark file name. Used later for subsetting pulled landmarks based on the chosen analysis type. 
get_spec_name <- function(landmark_path, position){
  clean_spec_name <- sub("_.{2}J_.+","", basename(landmark_path))
  return(clean_spec_name)
}

# Fit splines function
# apply_splines <- function(ridge_data) {
#   curve_list <- lapply(ridge_data, function(x) rbind(x, x[1, ]))
#   spline_applied <- lapply(curve_list, function(curve_ind) {
#     apply(curve_ind, 2, spline, method = "natural")
#   })
#   splines_binded <- lapply(spline_applied, function(spline_ind) {
#     y_values <- lapply(spline_ind, function(v) v$y)
#     do.call(cbind, y_values)
#   })
#   return(splines_binded)
# }

apply_splines <- function(ridge_data) {
  curve_list <- lapply(ridge_data, function(x) rbind(x, x[1, ]))
  # curve_list <- lapply(landmark_data_mats$EDJ_RIDGE, function(x) rbind(x, x[1, ]))
  
  spline_applied <- lapply(seq_along(curve_list), function(i) {
    curve_ind <- curve_list[[i]]
    # curve_ind <- curve_list[[1]]
    
    # Debug: Check for issues in this curve
    if(any(is.na(curve_ind)) | any(is.infinite(curve_ind))){
      print(paste("Processing curve", i))
      print(paste("Dimensions:", nrow(curve_ind), "x", ncol(curve_ind)))
      print(paste("Any NAs:", any(is.na(curve_ind))))
      print(paste("Any Inf:", any(is.infinite(curve_ind))))
      print(paste("WARNING: Column contains NA or Inf values. Exiting function."))
      return()
    }
    
    result <- apply(curve_ind, 2, function(col) {
      tryCatch(
        {
          spline(col, method = "natural")
        },
        warning = function(w) {
          message("Warning: ", conditionMessage(w))
          print(i)
          # still return something
          spline(col, method = "natural")
        },
        error = function(e) {
          message("Error: ", conditionMessage(e))
          NULL
        }
      )
    })
    
    
    return(result)
  })
  
  splines_binded <- lapply(spline_applied, function(spline_ind) {
    y_values <- lapply(spline_ind, function(v) v$y)
    do.call(cbind, y_values)
  })
  
  return(splines_binded)
}

# Get main landmark position for splitting the splines
get_main_landmark_position <- function(EDJ_spline, projected_main, iterant){
  d <- as.matrix(dist(rbind(projected_main[iterant, ], EDJ_spline)))
  fdistm <- matrix(d[1, -1]) 
  row.names(fdistm)<-c(1:length(fdistm))
  closestf<-as.numeric(names((fdistm[order(fdistm[,1]),])[c(1,2)]))
  
  if(iterant == 1){
    if(abs(diff(closestf))==1){
      fpos<-min(closestf)
    } else if(abs(diff(closestf))>1){
      fpos<- 0
    } 
  } else {
    fpos<-min(as.numeric(names((fdistm[order(fdistm[,1]),])[c(1,2)])))
  }
}

# Get plot boundaries for plotting from splines
spline_plot_boundaries <- function(splines){
  xmed<-median(splines[[a]][,1])
  ymed<-median(splines[[a]][,2])
  zmed<-median(splines[[a]][,3])
  return(list(xmed = xmed, ymed = ymed, zmed = zmed))
}

# Combine equidistant landmarks for each specimen
combine_equidistant_landmarks <- function(EDJ_eq = NULL, CEJ_eq = NULL, analysis){
  if(analysis == "CEJ"){
    combined <- eq_CEJ
  } else {
    nd <- length(EDJ_eq[[1]])  
    
    combined <- lapply(seq_len(nd), function(d) {
      EDJ_parts <- lapply(EDJ_eq, function(group) {
        mat <- group[[d]]
        mat[-nrow(mat), , drop = FALSE]  # remove last row
      })
      if(analysis == "EDJ CEJ"){
        CEJ_part <- eq_CEJ[[d]]
        do.call(rbind, c(EDJ_parts, list(CEJ_part)))
      } else {
        do.call(rbind, EDJ_parts)
      }
    })
  }
}

if(analysis == "EDJ CEJ"){
  file_pattern_o <- c(".EDJ_RIDGE.*.landmarkAscii", ".EDJ_MAIN.*.landmarkAscii", ".CEJ_RIDGE.*.landmarkAscii")
  # file_pattern <- paste(file_pattern, collapse ="|")
} else if(analysis == "EDJ"){
  file_pattern_o <- c(".EDJ_RIDGE.*.landmarkAscii", ".EDJ_MAIN.*.landmarkAscii")
  # file_pattern <- paste(file_pattern, collapse ="|")
} else if(analysis == "CEJ"){
  file_pattern_o <- c(".CEJ_RIDGE.*.landmarkAscii")
  # file_pattern <- paste(file_pattern, collapse ="|")
}

# Find all the landmarks on the server based on the provided position(s) and analysis type

landmark_files_list <- vector("list", length(file_pattern_o))
landmark_cats <- str_extract(file_pattern_o, ".{2}J_[^\\.]+")
names(landmark_files_list) <- landmark_cats

for(n in 1:length(landmark_files_list)){
  landmark_files_list[[n]] <- lapply(position, function(x) list.files(path = file.path("Landmarks", x), pattern = file_pattern_o[n], full.names = TRUE))
  names(landmark_files_list[[n]]) <- position
}

# Keep landmarks where the specimen name repeats across all landmark sets. 

# Extract the tooth type lists
for(pos in position){
  # Extract each sublist's Tooth type files
  um_lists <- lapply(landmark_files_list, `[[`, pos)
  
  # Extract IDs from each list
  id_lists <- lapply(um_lists, get_spec_name)
  
  # Get common IDs
  common_ids <- Reduce(intersect, id_lists)
  
  # Filter each sublist’s filenames
  filtered_lists <- mapply(function(files, ids) {
    files[get_spec_name(files) %in% common_ids]
  }, um_lists, id_lists, SIMPLIFY = FALSE)
  
  # Reassign filtered results back to file_list
  for (i in seq_along(landmark_files_list)) {
    landmark_files_list[[i]][[pos]] <- filtered_lists[[i]]
  }
}

# Collapse the tooth type lists
landmark_files_list_collapsed <- lapply(landmark_files_list, unlist, use.names = FALSE)

# Sort the landmarks
landmark_files_list_sorted <- lapply(landmark_files_list_collapsed, mixedsort)

# Get a single list of landmarks that will be used to construct a labels file
if(analysis  %in% c("CEJ", "EDJ CEJ")){
  spec_names <- get_spec_name(landmark_files_list_sorted$CEJ_RIDGE)
} else {
  spec_names <- get_spec_name(landmark_files_list_sorted$EDJ_RIDGE)
}

# load landmark files
landmark_data_lists <- lapply(landmark_files_list_sorted, function(file_group) {
  lapply(file_group, read.table, fill = TRUE, skip = 14, sep = "")
})

landmark_data_mats <- lapply(landmark_data_lists, function(land_table) {
  lapply(land_table, as.matrix)
})

# Apply splines based on analysis type
if (analysis %in% c("EDJ", "EDJ CEJ")) {
  EDJ_splines_binded <- apply_splines(landmark_data_mats$EDJ_RIDGE)
  
  for(c in 1:length(EDJ_splines_binded)){
    if(any(is.na(EDJ_splines_binded[[c]]))){
      print(c)
    }
  }
}

if (analysis %in% c("CEJ", "EDJ CEJ")) {
  CEJ_splines_binded <- apply_splines(landmark_data_mats$CEJ_RIDGE)
  
  for(c in 1:length(CEJ_splines_binded)){
    if(any(is.na(CEJ_splines_binded[[c]]))){
      print(c)
    }
  }
}

# Visualize the splines
a <- 1
clear3d("shapes")
if(analysis %in% c("EDJ CEJ", "EDJ")){
  plot_boundaries <- spline_plot_boundaries(EDJ_splines_binded)
  splines_to_plot <- EDJ_splines_binded
} else {
  plot_boundaries <- spline_plot_boundaries(CEJ_splines_binded)
  splines_to_plot <- CEJ_splines_binded
}

plot3d(splines_to_plot[[a]],type="l",xlim=c(plot_boundaries$xmed-10,plot_boundaries$xmed+10),ylim=c(plot_boundaries$ymed-10,plot_boundaries$ymed+10),zlim=c(plot_boundaries$zmed-10,plot_boundaries$zmed+10),aspect =T,xlab="",ylab="",zlab="",box=F,axes=F,main="Splines shown in black")

if(analysis %in% c("EDJ CEJ", "EDJ")){
  spheres3d(landmark_data_mats$CEJ_RIDGE[[a]],radius=0.05,col="blue")
  spheres3d(landmark_data_mats$EDJ_MAIN[[a]],radius=0.08,col="red")
  spheres3d(landmark_data_mats$EDJ_RIDGE[[a]],radius=0.05,col="blue")
  if(analysis == "EDJ CEJ"){
    lines3d(CEJ_splines_binded[[a]])
    texts3d(landmark_data_mats$CEJ_RIDGE[[a]],texts = c(0:length(landmark_data_mats$CEJ_RIDGE[[a]][,1])),font=2)
  }
  texts3d(landmark_data_mats$EDJ_RIDGE[[a]],texts = c(0:length(landmark_data_mats$EDJ_RIDGE[[a]][,1])),font=2)
} else {
  lines3d(CEJ_splines_binded[[a]])
  spheres3d(landmark_data_mats$CEJ_RIDGE[[a]],radius=0.08,col="red")
  texts3d(landmark_data_mats$CEJ_RIDGE[[a]],texts = c(0:length(landmark_data_mats$CEJ_RIDGE[[a]][,1])),font=2)
}


# Project main landmarks on to curve
if(analysis %in% c("EDJ CEJ", "EDJ")){
   proj_pts<-list()
  
  for(b in 1:length(landmark_data_mats$EDJ_RIDGE)){
    proj <- project_to_curve(landmark_data_mats$EDJ_MAIN[[b]],EDJ_splines_binded[[b]],stretch=0)
    proj_pts[[b]] <- proj$s
  }
  
  fposes_list <- vector("list", length = nrow(proj_pts[[1]]))
  
  fposes <- lapply(seq_along(fposes_list), function(l){
    mapply(function(EDJ_spline, projected_main){
      get_main_landmark_position(EDJ_spline, projected_main, iterant = l)
    }, 
    EDJ_splines_binded, proj_pts, SIMPLIFY = FALSE)
  })
  
  if(any(grepl("M", position))){
    number_of_f_points <- 4
    spline_list_name <- paste0("spline_", 1:4)
    spline_m_list_name <- paste0("spline_", 1:4,"m")
  } else {
    number_of_f_points <- 2
    spline_list_name <- paste0("spline_", 1:2)
    spline_m_list_name <- paste0("spline_", 1:2,"m")
  }
  
  # Initialize each f as an empty list
  EDJ_spline_segments <- vector("list", length = number_of_f_points)
  EDJ_spline_segments <- lapply(EDJ_spline_segments, function(x) list())
  names(EDJ_spline_segments) <- spline_list_name
  
  # For later placing the main landmarks at the end of each spline
  EDJ_spline_mains <- vector("list", length = number_of_f_points)
  EDJ_spline_mains <- lapply(EDJ_spline_mains, function(x) list())
  names(EDJ_spline_mains) <- spline_m_list_name
  
  for (i in seq_along(EDJ_splines_binded)) {
    f_vec <- sapply(fposes, "[[", i)
    spline <- EDJ_splines_binded[[i]]
    n <- nrow(spline)
    
    if (f_vec[[1]] < f_vec[[2]]) {
      lapply(seq_along(f_vec), function(f) {
        seg <- NULL
        if (length(f_vec) == 2) {
          if (f == 1) seg <- spline[(f_vec[[f]] + 1):f_vec[[f + 1]], ]
          else        seg <- spline[(f_vec[[f]] + 1):(n - 1), ]
        } else if (length(f_vec) == 4) {
          if (f == 1)       seg <- spline[(f_vec[[f]] + 1):f_vec[[f + 1]], ]
          else if (f < 4)   seg <- spline[(f_vec[[f]] + 1):f_vec[[f + 1]], ]
          else if (f == 4)  seg <- spline[(f_vec[[f]] + 1):(n - 1), ]
        }
        EDJ_spline_segments[[f]][[i]] <<- seg  # append matrix as ith element
      })
    } else {
      lapply(seq_along(f_vec), function(f) {
        seg <- NULL
        if (length(f_vec) == 2) {
          if (f == 1) seg <- spline[1:f_vec[[f + 1]], ]
          else        seg <- spline[c(f_vec[[f]] + 1):(f_vec[[f - 1]]), ]
        } else if (length(f_vec) == 4) {
          if (f == 1)       seg <- spline[1:f_vec[[f + 1]], ]
          else if (f < 4)   seg <- spline[(f_vec[[f]] + 1):f_vec[[f + 1]], ]
          else if (f == 4)  seg <- spline[(f_vec[[f]] + 1):f_vec[[f - 3]], ]
        }
        EDJ_spline_segments[[f]][[i]] <<- seg  # append matrix as ith element
      })
    }
  }
  
  # places main landmarks at end of each EDJ spline
  for (s in seq_along(EDJ_spline_mains)) {
    first_idx <- s
    last_idx  <- if (s < length(EDJ_spline_mains)) s + 1 else 1
    
    for (e in seq_along(landmark_data_mats$EDJ_MAIN)) {
      EDJ_spline_mains[[s]][[e]] <- rbind(
        landmark_data_mats$EDJ_MAIN[[e]][first_idx, ],
        EDJ_spline_segments[[s]][[e]],
        landmark_data_mats$EDJ_MAIN[[e]][last_idx, ]
      )
    }
  }
  
  a<-3
  {
    clear3d("shapes");

    if(analysis == "EDJ CEJ"){
      plot_boundaries <- spline_plot_boundaries(CEJ_splines_binded)
      splines_to_plot <- CEJ_splines_binded
    } else {
      plot_boundaries <- spline_plot_boundaries(EDJ_splines_binded)
      splines_to_plot <- EDJ_splines_binded
    }

    plot3d(splines_to_plot[[a]],type="l",xlim=c(plot_boundaries$xmed-10,plot_boundaries$xmed+10),ylim=c(plot_boundaries$ymed-10,plot_boundaries$ymed+10),zlim=c(plot_boundaries$zmed-10,plot_boundaries$zmed+10),aspect =T,xlab="",ylab="",zlab="",box=F,axes=F)
    if(any(grepl("P", position))){
      lines3d(EDJ_spline_mains$spline_1m[[a]],col="orange")
      texts3d(c(EDJ_spline_segments$spline_1[[a]][as.integer(length(EDJ_spline_segments$spline_1[[a]][,1])/2),]),texts = "EDJ 1",font=2)
      lines3d(EDJ_spline_mains$spline_2m[[a]],col="red")
      texts3d(c(EDJ_spline_segments$spline_2[[a]][as.integer(length(EDJ_spline_segments$spline_2[[a]][,1])/2),]),texts = "EDJ 2",font=2)
    } else if(any(grepl("M", position))){
      lines3d(EDJ_spline_mains$spline_1m[[a]],col="orange")
      texts3d(c(EDJ_spline_segments$spline_1[[a]][as.integer(length(EDJ_spline_segments$spline_1[[a]][,1])/2),]),texts = "EDJ 1",font=2)
      lines3d(EDJ_spline_mains$spline_2m[[a]],col="red")
      texts3d(c(EDJ_spline_segments$spline_2[[a]][as.integer(length(EDJ_spline_segments$spline_2[[a]][,1])/2),]),texts = "EDJ 2",font=2)     
      lines3d(EDJ_spline_mains$spline_3m[[a]],col="blue")
      texts3d(c(EDJ_spline_segments$spline_3[[a]][as.integer(length(EDJ_spline_segments$spline_3[[a]][,1])/2),]),texts = "EDJ 3",font=2)
      lines3d(EDJ_spline_mains$spline_4m[[a]],col="green")
      texts3d(c(EDJ_spline_segments$spline_4[[a]][as.integer(length(EDJ_spline_segments$spline_4[[a]][,1])/2),]),texts = "EDJ 4",font=2)
    }
    # spheres3d(landmark_data_mats$EDJ_MAIN[[a]],radius=0.05)
    if(analysis == "EDJ CEJ"){
      lines3d(CEJ_splines_binded[[a]],col="black")
    }
  }
}


# calculate maximum number of semilandmarks possible based on number of points splines are defined by
if(analysis %in% c("EDJ CEJ", "CEJ")){
  spline_to_iterate <- CEJ_splines_binded 
  minCEJ <- vector()
} else {
  spline_to_iterate <- EDJ_spline_mains$spline_1m 
}

for(i in 1:length(spline_to_iterate)){
  if(analysis %in% c("EDJ CEJ", "EDJ")){
    minEDJ <- lapply(EDJ_spline_mains, function(spline_num) {
      lapply(spline_num, function(spline) length(spline[,1]))
    })
    if(analysis == "EDJ CEJ")
    minCEJ[i]<-length(CEJ_splines_binded[[i]][,1])
  } else {
    minCEJ[i]<-length(CEJ_splines_binded[[i]][,1])
  }
}


if(analysis %in% c("EDJ CEJ", "EDJ")){
  print("Minimum number of landmarks:")
  for(l in 1:length(minEDJ)){
    print(paste0(paste0("EDJ",l), ": ", min(unlist(minEDJ[[l]]))))
  }

  EDJ_nPoints <- vector("list", length = number_of_f_points)
  names(EDJ_nPoints) <- paste0("EDJ", 1:number_of_f_points)
  
  # define the number of points for each EDJ curve. EDJ1, EDJ2 etc...
  points_per_curve <- c(25, 
                        35,
                        30,
                        20)
  
  for(i in 1:number_of_f_points){
    EDJ_nPoints[[i]] <- points_per_curve[i]
  }
  
  # define the numer of points for the CEJ 
  if(analysis == "EDJ CEJ"){
    print(paste0("CEJ: ", min(minCEJ)))
    CEJ<-60
  }
} else {
  print(paste0("CEJ: ", min(minCEJ)))
  CEJ<-40
}

# places equidistant semi-landmarks along the splines
# combined<-list()

eq_CEJ <- NULL
eq_EDJ <- NULL

if(analysis %in% c("EDJ CEJ", "EDJ")){
  one_spline <- EDJ_spline_mains$spline_1m[[1]]
  npoints <- EDJ_nPoints$EDJ1[[1]]
  eq_EDJ <- mapply(function(spline_group, npoints_group) {
    mapply(function(one_spline, npoints) {
      digit.curves(start = one_spline[1,], curve = one_spline, nPoints = npoints, closed = FALSE)
    }, spline_group, npoints_group, SIMPLIFY = FALSE)
  }, EDJ_spline_mains, EDJ_nPoints, SIMPLIFY = FALSE)
  
  if(analysis == "EDJ CEJ"){
    for(d in 1:length(CEJ_splines_binded)){
      eq_CEJ[[d]]<-digit.curves(start=CEJ_splines_binded[[d]][1,],curve=CEJ_splines_binded[[d]],nPoints=CEJ,closed=T)
    }
  }
} else if(analysis == "CEJ"){
  for(d in 1:length(CEJ_splines_binded)){
    eq_CEJ[[d]]<-digit.curves(start=CEJ_splines_binded[[d]][1,],curve=CEJ_splines_binded[[d]],nPoints=CEJ,closed=T)
    # combined[[d]]<-rbind(eq_CEJ[[d]][1:(length(eq_CEJ[[d]][,1])),])
  }
}

# Write how many landmarks are used for each section
if(analysis == "EDJ CEJ"){
  combined_num_lands <- c(points_per_curve, CEJ)
} else if(analysis == "EDJ"){
  combined_num_lands <- points_per_curve
} else {
  combined_num_lands <- CEJ
}

num_landmarks_df <- data.frame(matrix(combined_num_lands, ncol = length(combined_num_lands)))

if(analysis %in% c("EDJ", "EDJ CEJ")){
  if(any(grepl("M", position))){
    num_landmarks_df_colnames <- paste0("EDJ", c(1:4))
  } else {
    num_landmarks_df_colnames <- paste0("EDJ", c(1:2))
  }
  if(analysis == "EDJ CEJ"){
    num_landmarks_df_colnames <- c(num_landmarks_df_colnames, "CEJ")
  }
} else {
  num_landmarks_df_colnames <- "CEJ"
}

colnames(num_landmarks_df) <-  num_landmarks_df_colnames
  
# combine landmarks
combined <- combine_equidistant_landmarks(eq_EDJ, eq_CEJ, analysis)
NewCoordsArray <- array(as.numeric(unlist(combined)), dim=c(nrow(combined[[1]]), ncol(combined[[1]]) ,length(combined)))

for(i in 1:dim(NewCoordsArray)[3]){
  if(any(is.na(NewCoordsArray[,,i]))){
    print(paste0("NA detected in the landmark coordinates of this specimens, check landmarks! ", spec_names[i]))
  }
}

a<-1
clear3d("shapes")
if(analysis %in% c("EDJ CEJ", "CEJ")){
  plot_boundaries <- spline_plot_boundaries(CEJ_splines_binded)
  splines_to_plot <- CEJ_splines_binded
} else {
  plot_boundaries <- spline_plot_boundaries(EDJ_splines_binded)
  splines_to_plot <- EDJ_splines_binded
}

plot3d(splines_to_plot[[a]],type="l",xlim=c(plot_boundaries$xmed-10,plot_boundaries$xmed+10),ylim=c(plot_boundaries$ymed-10,plot_boundaries$ymed+10),zlim=c(plot_boundaries$zmed-10,plot_boundaries$zmed+10),aspect =T,xlab="",ylab="",zlab="",box=F,axes=F,main="Splines shown in black")
if(analysis %in% c("EDJ CEJ", "EDJ")){
  lines3d(EDJ_spline_mains$spline_1m[[a]],col="orange")
  points3d(eq_EDJ$spline_1m[[a]])
  texts3d(c(EDJ_spline_mains$spline_1m[[a]][as.integer(length(EDJ_spline_mains$spline_1m[[a]][,1])/2),]),texts = "EDJ 1",font=2)
  lines3d(EDJ_spline_mains$spline_2m[[a]],col="red")
  points3d(eq_EDJ$spline_2m[[a]])
  texts3d(c(EDJ_spline_mains$spline_2m[[a]][as.integer(length(EDJ_spline_mains$spline_2m[[a]][,1])/2),]),texts = "EDJ 2",font=2)
  if(any(grepl("M", position))){
    lines3d(EDJ_spline_mains$spline_3m[[a]],col="blue")
    points3d(eq_EDJ$spline_3m[[a]])
    texts3d(c(EDJ_spline_mains$spline_3m[[a]][as.integer(length(EDJ_spline_mains$spline_3m[[a]][,1])/2),]),texts = "EDJ 3",font=2)
    lines3d(EDJ_spline_mains$spline_4m[[a]],col="green")
    points3d(eq_EDJ$spline_4m[[a]])
    texts3d(c(EDJ_spline_mains$spline_4m[[a]][as.integer(length(EDJ_spline_mains$spline_4m[[a]][,1])/2),]),texts = "EDJ 4",font=2)
  }
  spheres3d(landmark_data_mats$EDJ_MAIN[[a]],radius=0.05)
  if(analysis == "EDJ CEJ"){
    points3d(eq_CEJ[[a]])
  }
} else {
    plot3d(splines_to_plot[[a]],type="l",xlim=c(plot_boundaries$xmed-10,plot_boundaries$xmed+10),ylim=c(plot_boundaries$ymed-10,plot_boundaries$ymed+10),zlim=c(plot_boundaries$zmed-10,plot_boundaries$zmed+10),aspect =T,xlab="",ylab="",zlab="",box=F,axes=F,main="Splines shown in black")
    points3d(eq_CEJ[[a]])
}

# plot3d(NewCoordsArray[,,a],type="l",xlim=c(plot_boundaries$xmed-10,plot_boundaries$xmed+10),ylim=c(plot_boundaries$ymed-10,plot_boundaries$ymed+10),zlim=c(plot_boundaries$zmed-10,plot_boundaries$zmed+10),aspect =T,xlab="",ylab="",zlab="",box=F,axes=F,main="Splines shown in black")
# text3d(NewCoordsArray[,,a], texts = 1:dim(NewCoordsArray)[1],font=2,cex=1)

if(analysis %in% c("EDJ CEJ", "EDJ")){
  curves <- vector("list", length(EDJ_nPoints))
  npoint_vals <- unlist(EDJ_nPoints)
  npoint_vals_l <- length(npoint_vals)
  
  if(analysis == "EDJ CEJ"){
    fix <- c(1, sapply(1:(length(npoint_vals)), function(i) sum(npoint_vals[1:i]) + (i+1)))
  } else {
    fix <- c(1, sapply(1:(length(npoint_vals)-1), function(i) sum(npoint_vals[1:i]) + (i+1)))
  }
  
  slid<-c(1:length(NewCoordsArray[,1,1]))[-c(fix)]
  
  # offsets
  offsets <- if(npoint_vals_l == 2){
    rep(2, npoint_vals_l)
  } else {
    first_seq <- seq_along(npoint_vals[-1]) + 1
    c(first_seq,  first_seq[3])
  } 
  
  # cumulative ends
  ends <- cumsum(npoint_vals) + offsets
  
  # start indices
  starts <- c(1, head(ends, -1))
  
  # generate curves
  curves <- mapply(function(s, e) s:e, starts, ends, SIMPLIFY = FALSE)
  
  # special wrap-around rule
  if(npoint_vals_l == 2){
    curves[[2]] <- c(curves[[2]], 1)  # second curve wraps to 1
  } else {
    curves[[npoint_vals_l]] <- c(curves[[npoint_vals_l]], 1)  # last curve wraps to 1 in longer cases
  }
  
  names(curves) <- names(EDJ_nPoints)
  
  if(analysis == "EDJ CEJ"){
    curveCEJ <- c((sum(npoint_vals) + npoint_vals_l + 1):(sum(npoint_vals) + CEJ + npoint_vals_l + 1), (sum(npoint_vals) + npoint_vals_l + 1))
    curves[[npoint_vals_l + 1]] <- curveCEJ
    names(curves)[[npoint_vals_l + 1]] <- "CEJ1"
  }
} else {
  fix<-c(1)
  curveCEJ <- c(1:(CEJ+1),1)
  curves<-list(curveCEJ)
  names(curves)<-"CEJ1"
  slid<-c(1:length(NewCoordsArray[,1,1]))[-c(fix)]
} 

# landmark_guide<-list(fixed_landmarks=fix,semilandmarks=slid,curves=curves) # Can't find where it is used?

# EDJ_curves<-c(curves[[1]],curves[[2]],curves[[3]],curves[[4]]) # Can't find where it is used?

# Make a template with the correct specimen names
template <- data.frame(Name=spec_names,Tooth=NA,Class1=NA,Class2=NA,Classify=NA,Exclude=NA)
template$Tooth <- str_extract(template$Name, ".{2}$")
template$Classify <- "*"
template$Exclude <- "*"
template$Class1 <- "*"
template$Class2 <- str_extract(template$Name, "^[^_]+")
if(length(position) > 1){
  position_file <- paste(position, collapse = "_")
} else {
  position_file <- position
}

write.table(template,paste(position_file,"_labels_extra.csv",sep=""),col.names = T,row.names=F,sep =",",quote=F)

# write the number of landmarks used
write.table(num_landmarks_df, paste(position_file,"_n_landmarks.csv",sep=""),sep =" ", col.names = T, row.names = F, quote = F)

# FILL IN TEMPLATE WITH SPECIMEN INFO

# Load filled in template
# spec_num<-length(spec_names) # Can't find where it is used?
classifierraw <-read.csv(paste(position_file,"_labels.csv",sep=""),header=TRUE)
keep<-classifierraw$Exclude!="ex" 
classifier<-classifierraw[keep,]
classified<-classifier$Classify!="classify"
spec_names_kept <- spec_names[keep]

if(analysis %in% c("EDJ CEJ", "EDJ")){
  # is it EDJ_spline_segments and not mains that need to be used....?
  # EDJ_spline_mains - because we are projecting both the slide and main landmarks back, no?
  
  EDJ_spline_subs <- lapply(seq_along(EDJ_spline_mains), function(spline_group){
    sub <- EDJ_spline_mains[[spline_group]][keep]
    names(sub) <- spec_names[keep]
    sub
  })
  
  names(EDJ_spline_subs) <- paste0("EDJ", seq_along(EDJ_spline_subs))
  
  datamat_EDJ_msub <- landmark_data_mats$EDJ_MAIN[keep]

  if(analysis == "EDJ CEJ"){
    CEJ_spline_sub<-CEJ_splines_binded[keep]
    names(CEJ_spline_sub)<-spec_names[keep]
    splines <- c(EDJ_spline_subs, list(CEJ = CEJ_spline_sub))
  } else{
    splines <- EDJ_spline_subs
  }
} else {
  CEJ_spline_sub<-CEJ_splines_binded[keep]
  names(CEJ_spline_sub)<-spec_names[keep]
}

# Either "Tooth_position" or "Taxonomy"
analysis_type <- "Tooth_position"

if(analysis_type == "Tooth_position"){
  gp <- as.factor(paste(classifier$Tooth)) # Tooth type
  
  farbe <- c(
    "black",
             "#4daf4a", 
             "#e41a1c",
             "#1874CD",
              "#984ea3", 
             "orange",
             "#53868B",
             "#8B864E",
             "#EEE8CD",
             "darkslategray",
             "#EE6AA7",
             "#54FF9F"
             )
} else {
  gp <- as.factor(paste(classifier$Class1)) # Taxonomy
  
  # set specific colours to group
  group_colours <- data.frame(
    matrix(
      c("Aus", "#e41a1c",
        "STS52", "black",
        "UW105", "black",
        "AL199", "black",
        "Early Homo", "#4daf4a",
        "East Paranthropus", "#984ea3",
        "Pboi", "#984ea3",
        "Hnea", "#984ea3",
        "Prob", "#1874CD",
        "ProbDNH", "#1C86EE",
        "Paet", "#53868B",
        "Hnal","#8B864E",
        "Here", "darkslategray",
        "Ased", "#EE5C42",
        "Hsap", "orange"
      ),
      ncol = 2, byrow = T)
  )
  
  colnames(group_colours) <- c("group", "colour") 
  
  farbe <- c()
  for(i in 1:nlevels(gp)){
    gp_ind <- which(group_colours$group == levels(gp)[i])
    farbe <- c(farbe, group_colours$colour[gp_ind])
  }
}

dimnames(NewCoordsArray)[[3]] <- classifierraw$Name
SubCoordsArray<-NewCoordsArray[,,keep]

for(i in 1:dim(SubCoordsArray)[3]){
  if(any(is.na(SubCoordsArray[,,i]))){
    print(paste0("NA detected in the landmark coordinates of this specimens, check landmarks! ", spec_names_kept[i]))
  }
}

# Find outliers - optional
# outliers <- find.outliers(SubCoordsArray)

### dropping landmarks 
# Plot landmarks of a specimen
a<-3;clear3d("shapes"); xmed<-median(SubCoordsArray[,,a][,1]); ymed<-median(SubCoordsArray[,,a][,2]); zmed<-median(SubCoordsArray[,,a][,3])
plot3d(SubCoordsArray[curves[[3]],,a],xlim=c(xmed-10,xmed+10),ylim=c(ymed-10,ymed+10),zlim=c(zmed-10,zmed+10),aspect =T,col="gray",type="l",xlab="",ylab="",zlab="",box=F,axes=F)
lines3d(SubCoordsArray[c(curves[[1]],curves[[2]]),,a],col="gray")
texts3d(SubCoordsArray[,,a],texts = c(1:length(SubCoordsArray[,,a][,1])),font=2)

# Plot PLY or STL with landmarks 
SurfaceModel <- file2mesh("O:/Research/DentalTissues/FOSSIL/NATIONAL_MUSEUMS_KENYA/KNM-ER820/KNM-ER820_10000258_002_001/KNM-ER 820_10000258_002_001_MAND.stl")
# SurfaceModel <- readSTL("O:/Research/DentalTissues/FOSSIL/NATIONAL_MUSEUMS_KENYA/KNM-ER820/KNM-ER820_10000258_002_001/KNM-ER 820_10000258_002_001_MAND.stl")
# SurfaceModel <- vcgPlyRead("O:/Research/DentalTissues/FOSSIL/UW_105/X0003973_000_001_UW_105-825_LRM/FILT1_1/UW105-825_LRM2-files/UW105-825_LLM2_Dentine.ply", updateNormals = TRUE, clean = TRUE)
SurfaceModel$material$color <- "Azure4"

dimnames(SubCoordsArray)
spec_index <- 48

plot3d(SurfaceModel, box = FALSE, aspect = FALSE, axes=FALSE, xlab="", ylab="", zlab="")
spheres3d(SubCoordsArray[,,spec_index], radius = 0.05, col = "red")
text3d(SubCoordsArray[,,spec_index], texts = 1:nrow(SubCoordsArray[,,spec_index]), col = "black", cex = 1.5, pos = 2)

# Drop landmarks
DROPLMs<-F
if (DROPLMs==TRUE) {
  #todrop<- c(1:4,19:26,41:47)
  todrop<-c(60:62, 1:6)
  #todrop<-c(1:46)
  total<-c(1:length(SubCoordsArray[,,a][,1]))
  tokeep<-total[-todrop]
} else {
  tokeep<-c(1:length(SubCoordsArray[,,a][,1]))
  todrop<-0
}

# First sliding step
slid1<-slider3d(dat.array=SubCoordsArray,SMvector=slid,outlines=curves,iterations = 1, stepsize = 1)

slid1_p<-array(dim=dim(slid1$dataslide))
for(a in 1:dim(slid1$dataslide)[3]){
  if(analysis %in% c("EDJ CEJ", "EDJ")){
    
    # Define segment boundaries and spline names dynamically
    if(analysis == "EDJ CEJ"){
      boundaries <- c(1, fix[2:length(fix)], dim(slid1$dataslide)[1] + 1)
      spline_refs <- c(paste0("EDJ", 1:(length(fix)-1)), "CEJ")
    } else {
      boundaries <- c(1, fix[2:length(fix)], dim(slid1$dataslide)[1] + 1)
      spline_refs <- paste0("EDJ", 1:(length(boundaries)-1))
    }
    
    # Use lapply to process each segment
    projected_list <- lapply(seq_along(spline_refs), function(i){
      seg_start <- boundaries[i]
      seg_end <- boundaries[i+1] - 1
      if(i == length(spline_refs)) seg_end <- dim(slid1$dataslide)[1] # Last segment
      
      segment <- seg_start:seg_end
      spline_name <- spline_refs[i]
      
      if(spline_name == "CEJ"){
        current_spline <- splines$CEJ[[a]]
      } else {
        current_spline <- splines[[spline_name]][[a]]
      }
      
      project_to_curve(slid1$dataslide[segment,,a], current_spline, stretch=0)$s
    })
    
    slid1_p[,,a] <- do.call(rbind, projected_list)
    slid1_p[fix,,a] <- slid1$dataslide[fix,,a]
    
  } else {
    s1_p_C <- project_to_curve(slid1$dataslide[c(fix[1]:dim(slid1$dataslide)[1]),,a], 
                               CEJ_spline_sub[[a]], stretch=0)
    slid1_p[,,a] <- s1_p_C$s
    slid1_p[fix,,a] <- slid1$dataslide[fix,,a]
  }
}

for(i in 1:dim(slid1_p)[3]){
  if(any(is.na(slid1_p[,,i]))){
    print(i)
  }
}

slid2<-slider3d(dat.array=slid1_p, SMvector=slid,outlines=curves,iterations = 1, stepsize = 1)

slid2_p<-array(dim=dim(slid1$dataslide))
for(a in 1:dim(slid1$dataslide)[3]){
  if(analysis %in% c("EDJ CEJ", "EDJ")){
    
    # Define segment boundaries and spline names dynamically
    if(analysis == "EDJ CEJ"){
      boundaries <- c(1, fix[2:length(fix)], dim(slid2$dataslide)[1] + 1)
      spline_refs <- c(paste0("EDJ", 1:(length(fix)-1)), "CEJ")
    } else {
      boundaries <- c(1, fix[2:length(fix)], dim(slid2$dataslide)[1] + 1)
      spline_refs <- paste0("EDJ", 1:(length(boundaries)-1))
    }
    
    # Use lapply to process each segment
    projected_list <- lapply(seq_along(spline_refs), function(i){
      seg_start <- boundaries[i]
      seg_end <- boundaries[i+1] - 1
      if(i == length(spline_refs)) seg_end <- dim(slid2$dataslide)[1] # Last segment
      
      segment <- seg_start:seg_end
      spline_name <- spline_refs[i]
      
      if(spline_name == "CEJ"){
        current_spline <- splines$CEJ[[a]]
      } else {
        current_spline <- splines[[spline_name]][[a]]
      }
      
      project_to_curve(slid2$dataslide[segment,,a], current_spline, stretch=0)$s
    })
    
    slid2_p[,,a] <- do.call(rbind, projected_list)
    slid2_p[fix,,a] <- slid2$dataslide[fix,,a]
    
  } else {
    s2_p_C <- project_to_curve(slid2$dataslide[c(fix[1]:dim(slid2$dataslide)[1]),,a], 
                               CEJ_spline_sub[[a]], stretch=0)
    slid2_p[,,a] <- s2_p_C$s
    slid2_p[fix,,a] <- slid2$dataslide[fix,,a]
  }
}

Landmarks_slid <- slid2_p

a<-3
clear3d("shapes"); 

if(analysis %in% c("EDJ CEJ", "EDJ")){
  # Detect available spline types dynamically
  edj_spline_names <- names(EDJ_spline_subs)[grepl("^EDJ", names(EDJ_spline_subs))]
  
  EDJ_spline_length <- length(EDJ_spline_subs[[edj_spline_names[1]]])
  
  # Function to get all splines for a given index
  get_splines_for_index <- function(i) {
    edj_splines <- setNames(lapply(edj_spline_names, function(name) EDJ_spline_subs[[name]][[i]]), 
                            edj_spline_names)
    edj_splines
  }
  
  # Get plot boundaries
  plot_boundaries <- if(analysis == "EDJ CEJ") {
    spline_plot_boundaries(CEJ_spline_sub)
  } else {
    spline_plot_boundaries(EDJ_spline_subs[[edj_spline_names[1]]])
  }
  
  # Create splines to plot
  splines_to_plot <- lapply(1:EDJ_spline_length, function(i){
    do.call(rbind, get_splines_for_index(i))
  })
  
  # Combined EDJ splines for lines3d
  combined_edj_splines <- lapply(1:EDJ_spline_length, function(i){
    do.call(rbind, lapply(edj_spline_names, function(name) EDJ_spline_subs[[name]][[i]]))
  })
  
} else {
  plot_boundaries <- spline_plot_boundaries(CEJ_spline_sub)
  splines_to_plot <- CEJ_spline_sub
}

# Rest of plotting code remains the same...
plot3d(splines_to_plot[[a]], type="l", xlab="", ylab="", zlab="", aspect=T, box=F, axes=F,
       xlim=c(plot_boundaries$xmed-10, plot_boundaries$xmed+10),
       ylim=c(plot_boundaries$ymed-10, plot_boundaries$ymed+10),
       zlim=c(plot_boundaries$zmed-10, plot_boundaries$zmed+10), col="gray")

if(analysis %in% c("EDJ CEJ", "EDJ")){
  lines3d(combined_edj_splines[[a]], col="gray")
  spheres3d(datamat_EDJ_msub[[a]], radius=0.1)
  if(analysis == "EDJ CEJ"){
    lines3d(CEJ_spline_sub[[a]], col="gray")
  }
} else {
  lines3d(CEJ_spline_sub[[a]], col="gray")
}

pch3d(SubCoordsArray[,,a], pch=19, cex=0.02, col="black")
pch3d(slid1$dataslide[,,a], pch=19, cex=0.02, col="#6767ff")
pch3d(slid1_p[,,a], pch=19, cex=0.04, col="#0000ff")
pch3d(slid2$dataslide[,,a], pch=19, cex=0.02, col="#ff7373")
pch3d(slid2_p[,,a], pch=19, cex=0.04, col="#ff0303")

n <- length(SubCoordsArray[,1,a])
Xslid1 <- rbind(SubCoordsArray[,,a], slid2_p[,,a])
OXslid1 <- Xslid1[as.vector(rbind(1:n, n + 1:n)), ]
try(arrows3d(OXslid1, scale=c(1,1,1), col="#ff0303"), silent=T)

Include_size<-F

# Procrustes registration
dimnames(Landmarks_slid)[3] <- dimnames(SubCoordsArray)[3]
Proc<-procSym(dataarray=Landmarks_slid,SMvector=NULL,outlines=curves, use.lm = tokeep, sizeshape = Include_size)

# Colours
# farbe<-c("black","#e41a1c","#984ea3","#377eb8","#4daf4a","black","black","#984ea3","black","black")

# Plot PCA PC1 V 2
plot(cbind(Proc$PCscores[,1],Proc$PCscores[,2]),type="n",asp=1,cex=1,
     xlab=paste0("PC 1 ", "(", round(Proc$Variance[1,2], 2), "%", ")"),
     ylab=paste0("PC 2 ", "(", round(Proc$Variance[2,2], 2), "%", ")"))
for(a in setdiff(seq_along(levels(gp)), c(0))){
  sub<-gp==levels(gp)[a]
  tr <- NULL
  try(tr<-tri.mesh(x=Proc$PCscores[sub,1],y=Proc$PCscores[sub,2],duplicate = "error"))
  if(!is.null(tr)){
    polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[a], alpha.f = 0.5)),border=farbe[a])
  } else if(sum(sub)==2){
    lines(Proc$PCscores[sub,1],Proc$PCscores[sub,2],col=farbe[a],lwd = 2)
  }}
points(cbind(Proc$PCscores[,1],Proc$PCscores[,2]),col=farbe[gp],pch=19)
# text(cbind(Proc$PCscores[,1],Proc$PCscores[,2]),label=dimnames(Proc$PCscores)[[1]],col=farbe[gp],pos=c(1,2),cex=0.6,offset=0.5)
# legend("bottomright", legend = levels(gp), col = farbe, pch = 19, bty = "n")

# Plot PCA PC1 V 3
plot(cbind(Proc$PCscores[,1],Proc$PCscores[,3]),type="n",asp=1,cex=1,
     xlab=paste0("PC 1 ", "(", round(Proc$Variance[1,2], 2), "%", ")"),
     ylab=paste0("PC 3 ", "(", round(Proc$Variance[3,2], 2), "%", ")"))
for(a in setdiff(seq_along(levels(gp)), c(0))){
  sub<-gp==levels(gp)[a]
  tr <- NULL
  try(tr<-tri.mesh(x=Proc$PCscores[sub,1],y=Proc$PCscores[sub,3],duplicate = "error"))
  if(!is.null(tr)){
    polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[a], alpha.f = 0.5)),border=farbe[a])
  } else if(sum(sub)==2){
    lines(Proc$PCscores[sub,1],Proc$PCscores[sub,3],col=farbe[a],lwd = 2)
  }}
points(cbind(Proc$PCscores[,1],Proc$PCscores[,3]),col=farbe[gp],pch=19)
text(cbind(Proc$PCscores[,1],Proc$PCscores[,3]),label=dimnames(Proc$PCscores)[[1]],col=farbe[gp],pos=c(3),cex=0.6,offset=0.5)
legend("topright", legend = levels(gp), col = farbe, pch = 19, bty = "n")

# 3D PCA
par3d(FOV = 0, windowRect = c(30, 30, 1770, 1000))

clear3d("shapes");plot.range<-1.02*max(abs(cbind(Proc$PCscores[,1],Proc$PCscores[,2],Proc$PCscores[,3])))
plot3d(cbind(Proc$PCscores[,1],Proc$PCscores[,2],Proc$PCscores[,3]),type="p",lwd=12,box=FALSE,size=10,col=farbe[gp],aspect =T,axes=T,
       xlim = c(min(Proc$PCscores[,1])+min(Proc$PCscores[,1])*0.10,max(Proc$PCscores[,1])+max(Proc$PCscores[,1])*0.10),
       ylim = c(min(Proc$PCscores[,2])+min(Proc$PCscores[,2])*0.10,max(Proc$PCscores[,2])+max(Proc$PCscores[,2])*0.10),
       zlim = c(min(Proc$PCscores[,3])+min(Proc$PCscores[,3])*0.10,max(Proc$PCscores[,3])+max(Proc$PCscores[,3])*0.10),
       xlab = paste("PC1 (", round(Proc$Variance[1,2],1), "%)", sep = ""),
       ylab = paste("PC2 (", round(Proc$Variance[2,2],1), "%)", sep = ""),
       zlab = paste("PC3 (", round(Proc$Variance[3,2],1), "%)", sep = ""))
# texts3d(cbind(Proc$PCscores[,1],Proc$PCscores[,2],Proc$PCscores[,3]),texts = rownames(Proc$PCscores),font=1,adj=1,cex=0.7)

for(a in setdiff(seq_along(levels(gp)), c(0))){
  sub<-gp==levels(gp)[a]
  PCsub<-cbind(Proc$PCscores[sub,1],Proc$PCscores[sub,2],Proc$PCscores[sub,3])
  if(length(PCsub[,1])>3){
    hull<-convhulln(PCsub)
    for(b in 1:length(hull[,1])){
      sub1<-PCsub[hull[b,],]
      triangles3d(sub1[,1],sub1[,2],sub1[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
    } 
  } else if(length(PCsub[,1])==3){
    triangles3d(PCsub[,1],PCsub[,2],PCsub[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
  } else if(length(PCsub[,1])==2){
    lines3d(PCsub,col=farbe[a])
  }
}
legend3d("topleft", levels(gp), pch = 16, col = farbe, cex = 3)

#Save 3D PCA plot as HTML file
save_location <- "D:/3D_GM/Projects/UW105/Results/Tooth_type_classifications/UW105_750_1007"
setwd(save_location)
widget <- rglwidget(webgl = TRUE, snapshot = FALSE, width = 1770, height = 1000)
htmlwidgets::saveWidget(widget, file.path(save_location, "Only_Paranthropus_clean.html"))

# Save 3D PCA as svg
rgl.postscript("Only_Paranthropus_clean.svg",fmt="svg")
setwd("D:/3D_GM/Projects/UW105")

# as.matrix(dist(Proc$rotated))

M <- mshape(Proc$rotated[tokeep,,])
PC <-Proc$PCscores[,1]
preds <- shape.predictor(Proc$rotated[tokeep,,], x= PC, Intercept = FALSE,
                         pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically

clear3d("shapes");
# xmed<-median(pos[,,PC][,1])
# ymed<-median(pos[,,PC][,2])
# zmed<-median(pos[,,PC][,3])

# min PC extreme
plot3d(rbind(preds$pred2, preds$pred2[1,]),type="l",xlab="", ylab="", zlab="",aspect =T,box=F,axes=F,xlim=c(xmed-0.4,xmed+0.4),ylim=c(ymed-0.4,ymed+0.4),zlim=c(zmed-0.4,zmed+0.4), lwd = 2)
spheres3d(preds$pred2[1,], radius=0.003)

clear3d("shapes");
# pos PC extreme
plot3d(rbind(preds$pred1, preds$pred1[1,]),type="l",xlab="", ylab="", zlab="",aspect =T,box=F,axes=F,xlim=c(xmed-0.4,xmed+0.4),ylim=c(ymed-0.4,ymed+0.4),zlim=c(zmed-0.4,zmed+0.4), lwd = 2)
spheres3d(preds$pred1[1,], radius=0.003)

lines3d(preds$pred1,col="red", lwd = 2)
text3d(preds$pred1, texts = 1:nrow(preds$pred1))


# Shape changes by PCs Morpho
pos<-pcaplot3d(Proc, pcshow = c(1,2,3), mag = 1, legend = F)
neg<-pcaplot3d(Proc, pcshow = c(1,2,3), mag = -1, legend = F)

PC<-1

clear3d("shapes");
xmed<-median(pos[,,PC][,1]); 
ymed<-median(pos[,,PC][,2]); 
zmed<-median(pos[,,PC][,3])
plot3d(pos[,,PC],type="l",xlab="", ylab="", zlab="",aspect =T,box=F,axes=F,xlim=c(xmed-0.4,xmed+0.4),ylim=c(ymed-0.4,ymed+0.4),zlim=c(zmed-0.4,zmed+0.4))
lines3d(neg[,,PC],col="red")

# nearest neighbour based on procrustes distance
dist_mat <- as.matrix(dist(two.d.array(Proc$rotated)))
rownames(dist_mat)
sort(dist_mat[which(rownames(dist_mat) == "UW105_966"),])

# compare specimens
as.matrix(classifier$Name)
#specimens<-c(1,17); cbind(as.character(classifier$Name[c(specimens)]),farbe[c(1:length(specimens))])
specimens<-c(1,16); cbind(as.character(classifier$Name[c(specimens)]),farbe[c(1:length(specimens))])
show_landmarks = F

clear3d("shapes");xmed<-median(Proc$rotated[,1,specimens[1]]); ymed<-median(Proc$rotated[,2,specimens[1]]); zmed<-median(Proc$rotated[,3,specimens[1]]);
plot3d(Proc$rotated[,,specimens[1]],type="n",xlab="", ylab="", zlab="",aspect =T,box=F,axes=F,xlim=c(xmed-0.4,xmed+0.4),ylim=c(ymed-0.4,ymed+0.4),zlim=c(zmed-0.4,zmed+0.4),col="red")

for(i in 1:length(specimens)){
  lapply(seq(curves), function(x){
    lines3d(Proc$rotated[curves[[x]],,specimens[i]], col=farbe[i],lwd=2)
    if(show_landmarks){
      points3d(Proc$rotated[curves[[x]],,specimens[i]],col=farbe[i], size = 7)
    }
    spheres3d(Proc$rotated[fix,,specimens[i]],radius=0.003, col=farbe[i])
  })
}
legend3d("topleft", cbind(as.character(classifier$Name[c(specimens)])), pch = 16, col = farbe, cex = 3)


## compare group means
# optional modification of the groups
# gp_char <- as.character(gp)
# gp_char[2] <- "UW105_966_mykolas"
# gp <- as.factor(gp_char)
# farbe <- c("black", farbe)

cbind(levels(gp),farbe[c(1:length(levels(gp)))])
groups<-c(1,3)

clear3d("shapes");
plot3d(Proc$rotated[,,1],type="n",xlab="", ylab="", zlab="",aspect =F,box=F,axes=F,col="red")
for(i in 1:length(groups)){
  sub<-gp==levels(gp)[groups[i]]
  group_mshape <- mshape(Proc$rotated[,,sub])
  
  lapply(seq(curves), function(x){
    lines3d(group_mshape[curves[[x]],], col=farbe[groups[i]],lwd=2)
    spheres3d(group_mshape[fix,],radius=0.003, col=farbe[groups[i]])
  })
}

legend3d("topleft", levels(gp)[groups], pch = 16, col = farbe[groups], cex = 3)

# Compare groups means + plots every specimens 
transparency_range_inv <- function(x, range_min, range_max){   
  (1 - (x - min(x))/(max(x)-min(x))) * (range_max - range_min) + range_min
}

nearest_neighbours_with_alpha <- function(procs, group_number, range_min = 0, range_max = 1){
  dist_mat <- as.matrix(dist(two.d.array(procs)))
  
  spec_match <- which(levels(gp)[group_number] == gp)
  spec_of_interest <- classifier$Name[spec_match]
  
  if(length(spec_of_interest) > 1){
    print("more than one specimen found to belong to the group!")
  } else if(length(spec_of_interest) == 0){
    print("no specimen was found to belong to the group!")
  }
  
  spec_idx <- which(rownames(dist_mat) == spec_of_interest)
  if (length(spec_idx) != 1)
    stop("Specimen name not found in distance matrix row names")
  
  
  dists <- dist_mat[spec_idx,]
  dists <- dists[-spec_idx]
  
  alpha_vals <- transparency_range_inv(dists, range_min , range_max)
  
  results_list <- list(distances = dists,
                       alpha_vals = alpha_vals)
}


dist_and_alpha_list<- nearest_neighbours_with_alpha(Proc$rotated, 1, 0.1 , 0.8)

gp_char <- as.character(gp)
gp_char[12] <- "UW105-1007_ULP4"
gp <- as.factor(gp_char)
farbe <- c(farbe, "black")

cbind(levels(gp),farbe[c(1:length(levels(gp)))])
groups<-c(1,2)

show_specimens = T

par3d(FOV = 0)
clear3d("shapes");plot3d(Proc$rotated[,,1],type="n",xlab="", ylab="", zlab="",aspect =F,box=F,axes=F,col="red")

for(i in 1:length(groups)){
  sub <- gp == levels(gp)[groups[i]]
  groups_specs <- Proc$rotated[,,sub]
  group_mshape <- mshape(groups_specs)
  
  lapply(seq(curves), function(x){
    lines3d(group_mshape[curves[[x]],], col=farbe[i], lwd=2)
    if(sum(sub) > 1 & show_specimens){
      for(s in 1:dim(groups_specs)[[3]]){
        lines3d(groups_specs[curves[[x]],,s], col=farbe[i], lwd=1.2, alpha=0.4)  
      }
    }
    spheres3d(group_mshape[fix,], radius=0.003, col=farbe[i])  
  })
}
legend3d("topleft", levels(gp)[groups], pch=16, col=farbe[1:length(groups)], cex=3)  

#Centroid size boxplot ggplot
library(ggrepel); library(ggplot2); library(viridis); library(tidyverse)
library(hrbrthemes)

gp <- factor(gsub("ProbDNH", "Prob",as.character(gp)))

manual_CEJ_csize <- apply(Landmarks_slid[unique(curves[length(curves)][[1]]),,], 3, cSize)
SizeFromProc <- Proc$size

# data.cs <- data.frame(CSize = log(Proc$size), Group = gp, Spec = rownames(Proc$PCscores))
data.cs <- data.frame(CSize = log(Proc$size), Group = gp, Spec = rownames(Proc$PCscores))

df_cs_box <- data.cs %>% filter(Spec != "BK_LM")
df_cs_line <- data.cs %>% filter(Spec == "BK_LM")


df_cs_box %>%
  ggplot(aes(x=Group, y=CSize, fill=Group)) +
  theme_bw() + 
  geom_boxplot(alpha = 0.5) +
  ggtitle("Log centroid size") +
  geom_point(size=2.5) + 
  geom_hline(data = df_cs_line, aes(yintercept = CSize),
             color = "black", linetype = "dashed", size = 0.6) +
  # geom_text_repel(aes(label = Spec), color = "black", size = 3, segment.color = "grey") +
  scale_color_manual(values=farbe) +
  scale_fill_manual(values=farbe) + 
  xlab("") + 
  ylab("") +
  theme(legend.position = "none") +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))


# Save the data as an R file
ProcrustesCoordinates <- Proc$rotated
RawLandmarks <- SubCoordsArray
FirstSliding <- slid1_p
SecondSliding <- slid2_p
save(RawLandmarks, ProcrustesCoordinates, FirstSliding, SecondSliding, classifier, file = "Results/Tooth_type_classifications/UW105_966/UW105_966_Paranthropus_anetrior_dentition_sample.RData")

#### CVA ####
Proc$Variance

dev.off()

pdf(paste("D:/3D_GM/Projects/UW105/Results/Tooth_type_classifications/UW105_966/PCs_CVA_labels_test.pdf",sep=""),height=15,width=15)
par(mfrow = c(3, 3))

# start PC = 80%, end PC = 95%
start_PC <- 3
end_PC <- 7

if(DROPLMs){
  coords_excluded<- Proc$PCscores_drop[!classified,1:end_PC]
} else {
  coords_excluded<- Proc$PCscores[!classified,1:end_PC]
}

gp_class <- factor(gp[classified])
groupnum_c <- nlevels(gp_class)
spec_names_classified <- spec_names_kept[classified]
spec_names_to_classify <- spec_names_kept[!classified]

coords_to_project<-NULL
if(sum(!classified)==1){
  coords_to_project<-t(as.matrix(coords_excluded))
  row.names(coords_to_project)<-classifier[!classified,"Name"]
} else if(sum(!classified)>1){
  coords_to_project<-as.matrix(coords_excluded)
}

res_df_classified <- data.frame(matrix("", ncol = (end_PC-start_PC)+3, nrow = length(spec_names_classified)+1))
colnames(res_df_classified) <- c("Class", paste0(start_PC:end_PC, " PCs"), "ClassAccuracy")
res_df_classified$Class <- c(as.character(gp_class), "-")
rownames(res_df_classified) <- c(spec_names_classified, "CVA accuracy")

res_df_to_classify <- data.frame(matrix("", ncol = (end_PC-start_PC)+3, nrow = length(spec_names_to_classify)))
colnames(res_df_to_classify) <- c("Class", paste0( start_PC:end_PC, " PCs"), "ClassResults")
rownames(res_df_to_classify) <- spec_names_to_classify

for(i in start_PC:end_PC){
  print(i)
  if(DROPLMs){
    CVA_res<-CVA(dataarray=Proc$PCscores_drop[classified,1:i],groups=gp_class,weighting = TRUE, cv = TRUE, prior=rep(1/groupnum_c,groupnum_c))
  } else {
    CVA_res<-CVA(dataarray=Proc$PCscores[classified,1:i],groups=gp_class,weighting = TRUE, cv = TRUE, prior=rep(1/groupnum_c,groupnum_c))
  }

  class_acc<-round(sum(CVA_res$class==CVA_res$groups)/length(CVA_res$groups)*100, 2)
  clas_res <- c(as.character(CVA_res$class), as.character(class_acc))
  res_df_classified[,which(names(res_df_classified) == paste0(i, " PCs"))] <- clas_res
  
  test<-classify(CVA_res,newdata=coords_to_project[,1:i, drop = FALSE])
  res_df_to_classify[,which(names(res_df_classified) == paste0(i, " PCs"))] <- as.character(test$class)
  
  
  # Plot CVA and project classified
  if(dim(coords_to_project)[1] < 2){
    swept <- coords_to_project[,1:i] -  CVA_res$Grandm
  } else {
    swept <- sweep(coords_to_project[,1:i], 2, CVA_res$Grandm)
  }
  # swept <- coords_to_project[1:i]- CVA_res$Grandm
  projected <- swept %*% CVA_res$CV
  rownames(projected) <- classifier$Name[!classified]
  
  if(ncol(CVA_res$CVscores) == 1){
    next
  }
  
  #CVA 1 v 2 (labels)
  yrange<-sum(abs(range(CVA_res$CVscores[,2])));ymin<-min(CVA_res$CVscores[,2]);ymax<-max(CVA_res$CVscores[,2])
  Xrange<-sum(abs(range(CVA_res$CVscores[,1])));Xmin<-min(CVA_res$CVscores[,1]);Xmax<-max(CVA_res$CVscores[,1])

  plot(cbind(CVA_res$CVscores[,1],CVA_res$CVscores[,2]),type="n",asp=1,cex=1,ylim=c(ymin-0.25*yrange,ymax+0.25*yrange),xlab=paste("CV 1 (",round(CVA_res$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA_res$Var[2,2]),"%)",sep=""))
  title(main = paste0(i ," PCs"))
  #plot(cbind(CVA$CVscores[,1],CVA$CVscores[,2]),type="n",asp=1,cex=1,xlim=c(Xmin-0.55*Xrange,Xmax+0.35*Xrange),xlab=paste("CV 1 (",round(CVA$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA$Var[2,2]),"%)",sep=""))
  for(gruppe in 1:length(levels(gp_class))){
    sub<-gp_class==levels(gp_class)[gruppe]
    tr <- NULL
    try(tr<-tri.mesh(x=CVA_res$CVscores[sub,1],y=CVA_res$CVscores[sub,2],duplicate = "error"))
    if(!is.null(tr)){
      polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[gruppe], alpha.f = 0.5)),border=farbe[gruppe])
    } else if(sum(sub)==2){
      lines(CVA_res$CVscores[sub,1],CVA_res$CVscores[sub,2],col=farbe[gruppe],lwd = 2)
    }}
  points(CVA_res$CVscores[,c(1,2)],col=farbe[gp_class],pch=19)
  # text(CVA_res$CVscores[,c(1,2)],label=classifier$Name[classified],col=farbe[classifier$Class2[classified]],pos=c(1,2),cex=0.6,offset=0.5)
  legend("bottomright", legend = levels(gp_class), col = farbe, pch = 19, bty = "n")
  points(cbind(projected[,1],projected[,2]),pch=19)
  text(cbind(projected[,1],projected[,2]), label = rownames(projected), cex = 0.6, pos = c(1,2))
}


res_df_classified$ClassAccuracy <- apply(res_df_classified, 1, function(x){
  correct_freq <- x[2:(length(x)-1)] == x[1]
  x[length(x)] <- round(((sum(correct_freq)/length(correct_freq))*100), 2)
})

res_df_to_classify$ClassResults <- apply(res_df_to_classify, 1, function(x){
  class_tab <- table(x)
  class_tab_percent <- round(class_tab[2:length(class_tab)]/sum(class_tab[2:length(class_tab)])*100,2)
  class_res <- paste0(names(class_tab_percent), ": ", class_tab_percent, "%", collapse = " ")
  class_res
})

per_spec_total_class <- mean(res_df_classified$ClassAccuracy[-length(res_df_classified$ClassAccuracy)])

res_df_classified$ClassAccuracy[nrow(res_df_classified)] <- per_spec_total_class

CVA_results <- list(
  ClassifiedResults = res_df_classified,
  ToClassifyResults = as.data.frame(t(res_df_to_classify))
)

CVA_results$ToClassifyResults

class_mat <- data.frame(matrix(c(gp_class, CVA_res$class), ncol =2))
rownames(class_mat) <- spec_names_kept[classified]


# Plot CVA and project classified
swept <- sweep(coords_to_project, 2, CVA_res$Grandm)
projected <- swept %*% CVA_res$CV

swep_all <- sweep(Proc$PCscores_drop[classified,1:i], 2, CVA_res$Grandm)
projected_all <- swep_all %*% CVA_res$CV

#CVA 1 v 2 (labels)
yrange<-sum(abs(range(CVA_res$CVscores[,2])));ymin<-min(CVA_res$CVscores[,2]);ymax<-max(CVA_res$CVscores[,2])
Xrange<-sum(abs(range(CVA_res$CVscores[,1])));Xmin<-min(CVA_res$CVscores[,1]);Xmax<-max(CVA_res$CVscores[,1])

#pdf(paste("../Figures/All PCAs/Mar23 - Nat_comm_R1_Sang/CVA/",position,"_",number,"PCs_CVA_labels.pdf",sep=""),height=7,width=width)
plot(cbind(CVA_res$CVscores[,1],CVA_res$CVscores[,2]),type="n",asp=1,cex=1,ylim=c(ymin-0.25*yrange,ymax+0.25*yrange),xlab=paste("CV 1 (",round(CVA_res$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA_res$Var[2,2]),"%)",sep=""))
#plot(cbind(CVA$CVscores[,1],CVA$CVscores[,2]),type="n",asp=1,cex=1,xlim=c(Xmin-0.55*Xrange,Xmax+0.35*Xrange),xlab=paste("CV 1 (",round(CVA$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA$Var[2,2]),"%)",sep=""))
for(gruppe in 1:length(levels(gp_class))){
  sub<-gp_class==levels(gp_class)[gruppe]
  tr <- NULL
  try(tr<-tri.mesh(x=CVA_res$CVscores[sub,1],y=CVA_res$CVscores[sub,2],duplicate = "error"))
  if(!is.null(tr)){
    polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[gruppe], alpha.f = 0.5)),border=farbe[gruppe])
  } else if(sum(sub)==2){
    lines(CVA_res$CVscores[sub,1],CVA_res$CVscores[sub,2],col=farbe[gruppe],lwd = 2)
  }}
points(CVA_res$CVscores[,c(1,2)],col=farbe[gp_class],pch=19)
text(CVA_res$CVscores[,c(1,2)],label=classifier$Name[classified],col=farbe[classifier$Class2[classified]],pos=c(1,2),cex=0.6,offset=0.5)
points(cbind(projected[,1],projected[,2]),pch=19)
text(cbind(projected[,1],projected[,2]), label = rownames(projected), cex = 0.6, pos = c(1,2))

# typicality probabilities
res_df_typicality_prob <- data.frame(matrix("", ncol = (end_PC-start_PC)+1, nrow = length(spec_names_to_classify)))
colnames(res_df_typicality_prob) <- paste0(start_PC:end_PC, " PCs")
rownames(res_df_typicality_prob) <- spec_names_to_classify

for(i in start_PC:end_PC){
  if(DROPLMs){
    CVA_res<-CVA(dataarray=Proc$PCscores_drop[classified,1:i],groups=gp_class,weighting = TRUE, cv = TRUE, prior=rep(1/groupnum_c,groupnum_c))
  } else {
    CVA_res<-CVA(dataarray=Proc$PCscores[classified,1:i],groups=gp_class,weighting = TRUE, cv = TRUE, prior=rep(1/groupnum_c,groupnum_c))
  }

  if(dim(coords_to_project)[1] < 2){
    swept <- coords_to_project[,1:i] -  CVA_res$Grandm
  } else {
    swept <- sweep(coords_to_project[,1:i], 2, CVA_res$Grandm)
  }
  projected <- swept %*% CVA_res$CV
  typ<-typprobClass(x=projected,data=CVA_res$CVscores,groups=gp_class,method="wilson",outlier=0.01)
  res_df_typicality_prob[,which(names(res_df_typicality_prob) == paste0(i, " PCs"))] <- as.character(typ$groupaffin)
}

typicality_probs_results <- data.frame(t(res_df_typicality_prob))
typicality_probs_results

write.csv(typicality_probs_results, "Results/UW105-811_LM1_typicality_probs.csv", col.names = T, row.names = F, quote = F)

D2 <- mahalanobis(projected[,], 
                  colMeans(CVA_res$CVscores[gp_class == "Prob",]), 
                  cov(CVA_res$CVscores[gp_class == "Prob",]))
D2

typ_class <- typprobClass(x=CVA_res$CVscores,groups=gp_class,method="c",outlier=0.1)

# 3D PCA
par3d(FOV = 0, windowRect = c(30, 30, 1770, 1000))

clear3d("shapes");plot.range<-1.02*max(abs(cbind(CVA_res$CVscores[,1],CVA_res$CVscores[,2],CVA_res$CVscores[,3])))
plot3d(cbind(CVA_res$CVscores[,1],CVA_res$CVscores[,2],CVA_res$CVscores[,3]),type="p",lwd=12,box=FALSE,size=10,col=farbe[gp],aspect =T,axes=T,
       xlim = c(min(CVA_res$CVscores[,1])+min(CVA_res$CVscores[,1])*0.10,max(CVA_res$CVscores[,1])+max(CVA_res$CVscores[,1])*0.10),
       ylim = c(min(CVA_res$CVscores[,2])+min(CVA_res$CVscores[,2])*0.10,max(CVA_res$CVscores[,2])+max(CVA_res$CVscores[,2])*0.10),
       zlim = c(min(CVA_res$CVscores[,3])+min(CVA_res$CVscores[,3])*0.10,max(CVA_res$CVscores[,3])+max(CVA_res$CVscores[,3])*0.10),
       xlab = paste("PC1 (", round(Proc$Variance[1,2],1), "%)", sep = ""),
       ylab = paste("PC2 (", round(Proc$Variance[2,2],1), "%)", sep = ""),
       zlab = paste("PC3 (", round(Proc$Variance[4,2],1), "%)", sep = ""))
texts3d(cbind(CVA_res$CVscores[,1],CVA_res$CVscores[,2],CVA_res$CVscores[,3]),texts = classifier$Name[classified],font=2,adj=1,cex=0.7)
legend3d("topleft", levels(gp_class), pch = 16, col = farbe, cex = 3)

for(a in 1:length(levels(gp_class))){
  sub<-gp_class==levels(gp_class)[a]
  PCsub<-cbind(CVA_res$CVscores[sub,1],CVA_res$CVscores[sub,2],CVA_res$CVscores[sub,3])
  if(length(PCsub[,1])>3){
    hull<-convhulln(PCsub)
    for(b in 1:length(hull[,1])){
      sub1<-PCsub[hull[b,],]
      triangles3d(sub1[,1],sub1[,2],sub1[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
    } 
  } else if(length(PCsub[,1])==3){
    triangles3d(PCsub[,1],PCsub[,2],PCsub[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
  } else if(length(PCsub[,1])==2){
    lines3d(PCsub,col=farbe[a])
  }
}

points3d(cbind(projected[,1],projected[,2],projected[,3]), size =10)
texts3d(cbind(projected[,1],projected[,2],projected[,3]),texts = rownames(projected),font=2,adj=1,cex=0.7)
legend3d("topleft", levels(gp), pch = 16, col = farbe, cex = 3)


# CVA accuracy by PCs used

if(DROPLMs){
  ncomp <- ncol(Proc$PCscores_drop)
} else{
  ncomp <- ncol(Proc$PCscores)
}

dat <- as.data.frame(cbind(Group = gp_class, Proc$PCscores[classified, 1:ncomp]))
accuracies <- rep(NA, ncomp - 1)
names(accuracies) <- paste(2:ncomp, "PCs")
for (k in 3:ncol(dat)) {
  
  #mod <- lda(dat[, 2:k],grouping=dat$Group,prior = rep(1/groupnum_c, groupnum_c),CV = TRUE)
  mod <- lda(dat$Group ~ ., data = dat[, 1:k], prior=rep(1/groupnum_c,groupnum_c), CV = TRUE)
  accuracies[k - 2] <- sum(mod$class == dat$Group,na.rm=T) / length(na.exclude(mod$class))
}

plot(x = 2:ncomp, y = accuracies,
     type = "b", pch = 16,
     xlab = "Number of PCA axes in LDA model",
     ylab = "X-val classification accuracy",
     ylim = (c(0,1)))
grid()


# compare group means
cbind(levels(gp),farbe[c(1:length(levels(gp)))])
groups<-c(1,3,5)
groups<-c(1,2)

clear3d(c("shapes","background"));plot3d(Proc$rotated[,,1],type="n",xlab="", ylab="", zlab="",aspect =F,box=F,axes=F)
for(i in 1:length(groups)){
  sub<-gp==levels(gp)[groups[i]]
  X<-Proc$rotated[,1,sub]
  Y<-Proc$rotated[,2,sub]
  Z<-Proc$rotated[,3,sub]
  if(is.vector(X)){
    lines3d(Proc$rotated[c(curves[[1]],curves[[2]],curves[[3]],curves[[4]]),,sub],col=farbe[groups[i]],lwd=2)
    #lines3d(Proc$rotated[curves[[3]],,sub],col=farbe[groups[i]],lwd=2)
    spheres3d(Proc$rotated[fix[c(1,2,3,4)],,sub],radius=0.002, col=farbe[groups[i]])
  } else {
    mean<-cbind(apply(X,1,mean),apply(Y,1,mean),apply(Z,1,mean))
    lines3d(mean[c(curves[[1]],curves[[2]],curves[[3]],curves[[4]]),],col=farbe[groups[i]],lwd=2)
    #lines3d(mean[curves[[3]],],col=farbe[groups[i]],lwd=2)
    spheres3d(mean[fix[c(1,2,3,4)],],radius=0.002, col=farbe[groups[i]])
  }
}

legend3d("topright", legend = levels(gp)[groups], pch = 16, col = farbe[groups], cex=1, inset=c(0.02))

##REGENERATE PCA ON LM SUBSET
include_size <- F

procrot <- Morpho:::orp(Proc$rotated[tokeep,,], mshape = Proc$mshape[tokeep,])
Symtan <- procrot
Symtan <- sweep(Symtan, 1:2, Proc$mshape[tokeep,])
tan <- vecx(Symtan)
if(include_size){
  csizes = log(apply(SubCoordsArray[tokeep,,],3,cSize))
  princ <- prcompfast(cbind(tan, csizes), scale. = FALSE)
} else {
  princ <- prcompfast(tan)
}

values <- 0
eigv <- princ$sdev^2
values <- eigv[which(eigv > 1e-14)]
lv <- length(values)
PCs <- princ$rotation[, 1:lv]
PCscore_sym <- as.matrix(princ$x[, 1:lv])
Proc$PCscores_drop<-PCscore_sym

if (length(values) == 1) {
  SymVar <- values
} else {
  SymVar <- matrix(NA, length(values), 3)
  SymVar[, 1] <- values
  for (i in 1:length(values)) SymVar[i, 2] <- (values[i]/sum(values)) * 
    100
  SymVar[1, 3] <- SymVar[1, 2]
  for (i in 2:length(values)) SymVar[i, 3] <- SymVar[i, 
                                                     2] + SymVar[i - 1, 3]
  colnames(SymVar) <- c("eigenvalues", "% Variance", "Cumulative %")
  rownames(SymVar) <- 1:nrow(SymVar)
}

Proc$Variance_drop<-SymVar

# Plot PCA PC1 V 2
plot(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,2]),type="n",asp=1,cex=1,
     xlab=paste0("PC 1 ", "(", round(Proc$Variance[1,2], 2), "%", ")"),
     ylab=paste0("PC 2 ", "(", round(Proc$Variance[2,2], 2), "%", ")"))
for(a in setdiff(seq_along(levels(gp)), c(6))){
  sub<-gp==levels(gp)[a]
  tr <- NULL
  try(tr<-tri.mesh(x=Proc$PCscores_drop[sub,1],y=Proc$PCscores_drop[sub,2],duplicate = "error"))
  if(!is.null(tr)){
    polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[a], alpha.f = 0.5)),border=farbe[a])
  } else if(sum(sub)==2){
    lines(Proc$PCscores_drop[sub,1],Proc$PCscores_drop[sub,2],col=farbe[a],lwd = 2)
  }}
points(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,2]),col=farbe[gp],pch=19)
text(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,2]),label=dimnames(SubCoordsArray)[[3]],col=farbe[gp],pos=c(1,2),cex=0.6,offset=0.5)
legend("topright", legend = levels(gp), col = farbe, pch = 19, bty = "n")

# Plot PCA PC1 V 3
plot(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,3]),type="n",asp=1,cex=1,
     xlab=paste0("PC 1 ", "(", round(Proc$Variance[1,2], 2), "%", ")"),
     ylab=paste0("PC 3 ", "(", round(Proc$Variance[3,2], 2), "%", ")"))
for(a in setdiff(seq_along(levels(gp)), c(0))){
  sub<-gp==levels(gp)[a]
  tr <- NULL
  try(tr<-tri.mesh(x=Proc$PCscores_drop[sub,1],y=Proc$PCscores_drop[sub,3],duplicate = "error"))
  if(!is.null(tr)){
    polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[a], alpha.f = 0.5)),border=farbe[a])
  } else if(sum(sub)==2){
    lines(Proc$PCscores_drop[sub,1],Proc$PCscores_drop[sub,3],col=farbe[a],lwd = 2)
  }}
points(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,3]),col=farbe[gp],pch=19)
# text(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,3]),label=dimnames(SubCoordsArray)[[3]],col=farbe[gp],pos=c(1,2),cex=0.6,offset=0.5)
# legend("topright", legend = levels(gp), col = farbe, pch = 19, bty = "n")

# 3D PCA

par3d(FOV = 0, windowRect = c(30, 30, 1770, 1000))

clear3d("shapes");plot.range<-1.02*max(abs(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,2],Proc$PCscores_drop[,3])))
plot3d(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,2],Proc$PCscores_drop[,3]),type="p",lwd=12,box=FALSE,size=10,col=farbe[gp],aspect =T,axes=T,
       xlim = c(min(Proc$PCscores_drop[,1])+min(Proc$PCscores_drop[,1])*0.10,max(Proc$PCscores_drop[,1])+max(Proc$PCscores_drop[,1])*0.10),
       ylim = c(min(Proc$PCscores_drop[,2])+min(Proc$PCscores_drop[,2])*0.10,max(Proc$PCscores_drop[,2])+max(Proc$PCscores_drop[,2])*0.10),
       zlim = c(min(Proc$PCscores_drop[,3])+min(Proc$PCscores_drop[,3])*0.10,max(Proc$PCscores_drop[,3])+max(Proc$PCscores_drop[,3])*0.10),
       xlab = paste("PC1 (", round(Proc$Variance[1,2],1), "%)", sep = ""),
       ylab = paste("PC2 (", round(Proc$Variance[2,2],1), "%)", sep = ""),
       zlab = paste("PC3 (", round(Proc$Variance[3,2],1), "%)", sep = ""))
# texts3d(cbind(Proc$PCscores_drop[,1],Proc$PCscores_drop[,2],Proc$PCscores_drop[,3]),texts = classifier$Name,font=2,adj=1,cex=0.7)
legend3d("topleft", levels(gp), pch = 16, col = farbe, cex = 3)

for(a in setdiff(seq_along(levels(gp)), c(6))){
  sub<-gp==levels(gp)[a]
  PCsub<-cbind(Proc$PCscores_drop[sub,1],Proc$PCscores_drop[sub,2],Proc$PCscores_drop[sub,3])
  if(length(PCsub[,1])>3){
    hull<-convhulln(PCsub)
    for(b in 1:length(hull[,1])){
      sub1<-PCsub[hull[b,],]
      triangles3d(sub1[,1],sub1[,2],sub1[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
    } 
  } else if(length(PCsub[,1])==3){
    triangles3d(PCsub[,1],PCsub[,2],PCsub[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
  } else if(length(PCsub[,1])==2){
    lines3d(PCsub,col=farbe[a])
  }
}
legend3d("topleft", levels(gp), pch = 16, col = farbe, cex = 3)

#Save 3D PCA plot as HTML file
save_location <- "D:/3D_GM/Projects/UW105/Results/Taxonomic_assessment/UW105_240_811_LP4"
setwd(save_location)
widget <- rglwidget(webgl = TRUE, snapshot = FALSE, width = 1770, height = 1000)
htmlwidgets::saveWidget(widget, file.path(save_location, "Shape_PCA_3D_dropped_landmarks_clean.html"))

# Save 3D PCA as svg
rgl.postscript("Shape_PCA_3D_dropped_landmarks_clean.svg",fmt="svg")
setwd("D:/3D_GM/Projects/UW105")

######### 
library(ggrepel); library(ggplot2); library(viridis); library(tidyverse)
library(hrbrthemes)

data.cs <- data.frame(CSize = apply(SubCoordsArray[tokeep,,],3,cSize), Group = gp, Spec = rownames(Proc$PCscores))

df_cs_box <- data.cs %>% filter(Spec != "CC7_2045")
df_cs_line <- data.cs %>% filter(Spec == "CC7_2045")

df_cs_box %>%
  ggplot(aes(x=Group, y=CSize, fill=Group)) +
  theme_bw() + 
  geom_boxplot(alpha = 0.5) +
  ggtitle("Log centroid size") +
  geom_point(size=2.5) + 
  geom_hline(data = df_cs_line, aes(yintercept = CSize),
             color = "black", linetype = "dashed", size = 0.6) +
  # geom_text_repel(aes(label = Spec), color = "black", size = 3, segment.color = "grey") +
  scale_color_manual(values=farbe[-1]) +
  scale_fill_manual(values=farbe[-1]) + 
  xlab("") + 
  ylab("") +
  theme(legend.position = "none") +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=15),
        plot.title = element_text(hjust = 0.5))

## premolar group comparisons 


# DROP LANDMARKS
DROPLMs<-TRUE
if (DROPLMs==TRUE) {
  
  # todrop<- c(63:67, 81:84, 1:6)
  
  total<-c(1:length(SubCoordsArray[,,1][,1]))
  tokeep<-total[-todrop]
  
  curves_drop_orgn <- lapply(seq(curves), function(x){
    kept_curve <- curves[[x]][!curves[[x]] %in% todrop]
  })
  
  ds <- lapply(seq(curves_drop_orgn), function(x){
    d <- length(curves_drop_orgn[[x]])
  })
  
  replacement<-cbind(c(1:length(tokeep)),tokeep)
  fix_drop<-replacement[match(fix,replacement[,2]),1]
  fix_drop<-fix_drop[!is.na(fix_drop)]
  curves_drop<-list()
  for(i in 1:length(curves_drop_orgn)){
    curves_drop[[i]]<-replacement[match(curves_drop_orgn[[i]],replacement[,2]),1]
  }
  warning("curves_drop will not work correctly if drops are in middle of a curve")
} else {
  tokeep<-c(1:length(SubCoordsArray[,,1][,1]))
  todrop<-0
  curves_drop_orgn<-curves;curves_drop<-curves;fix_drop<-fix
}

# compare between groups
include_size<-FALSE

par3d(FOV = 0)

cbind(levels(gp),farbe[c(1:length(levels(gp)))])
# groups<-c(1,3,5)
groups<-c(2,6)

clear3d("shapes");
plot3d(Proc$rotated[curves_drop_orgn[[1]],,1],type="n",xlab="", ylab="", zlab="",aspect =F,box=F,axes=F,col="red")
for(i in 1:length(groups)){
  sub<-gp==levels(gp)[groups[i]]
  group_mshape <- mshape(Proc$rotated[,,sub])
  
  lapply(seq(curves_drop_orgn), function(x){
    lines3d(group_mshape[curves_drop_orgn[[x]],], col=farbe[groups[i]],lwd=2)
    match_fixed <- curves_drop_orgn[[x]] %in% fix 
    if(any(match_fixed)){
      spheres3d(group_mshape[curves_drop_orgn[[x]][match_fixed],],radius=0.003, col=farbe[groups[i]])
    }
    # text3d(Proc$rotated[curves_drop_orgn[[x]],,sub], texts = 1:length(curves_drop_orgn[[x]])) 
  })
}

legend3d("topleft", levels(gp)[groups], pch = 16, col = farbe[groups], cex = 3)

# 
# groups<-c(1,6)
# 
# clear3d("shapes");
# plot3d(Proc$rotated[curves_drop_orgn[[1]],,1],type="n",xlab="", ylab="", zlab="",aspect =F,box=F,axes=F,col="red")
# for(i in 1:length(groups)){
#   sub<-gp==levels(gp)[groups[i]]
#   X<-Proc$rotated[,1,sub]
#   Y<-Proc$rotated[,2,sub]
#   Z<-Proc$rotated[,3,sub]
#   if(is.vector(X)){
#     if(analysis %in% c("EDJ CEJ", "EDJ")){
#       lines3d(Proc$rotated[curves_drop_orgn[[1]],,sub],col=farbe[groups[i]],lwd=2)
#       lines3d(Proc$rotated[curves_drop_orgn[[2]],,sub],col=farbe[groups[i]],lwd=2)
#       lines3d(Proc$rotated[curves[[3]],,sub],col=farbe[groups[i]],lwd=2)
#     } else{
#       # lines3d(Proc$rotated[curves_drop_orgn[[1]],,sub],col=farbe[groups[i]],lwd=2)
#       lines3d(Proc$rotated[curves_drop_orgn[[1]][1:36],,sub],col=farbe[groups[i]],lwd=2)
#       lines3d(Proc$rotated[curves_drop_orgn[[1]][37:length(curves_drop_orgn[[1]])],,sub],col=farbe[groups[i]],lwd=2)
#     }
#     spheres3d(Proc$rotated[fix_drop,,sub],radius=0.003, col=farbe[groups[i]])
#   } else {
#     mean<-cbind(apply(X,1,mean),apply(Y,1,mean),apply(Z,1,mean))
#     if(analysis %in% c("EDJ CEJ", "EDJ")){
#       lines3d(mean[curves_drop_orgn[[1]],],col=farbe[groups[i]],lwd=2)
#       lines3d(mean[curves_drop_orgn[[2]],],col=farbe[groups[i]],lwd=2)
#       lines3d(mean[curves[[3]],],col=farbe[groups[i]],lwd=2)
#     } else{
#       lines3d(mean[curves_drop_orgn[[1]][1:36],],col=farbe[groups[i]],lwd=2)
#       lines3d(mean[curves_drop_orgn[[1]][37:length(curves_drop_orgn[[1]])],],col=farbe[groups[i]],lwd=2)
#     }
#     # lines3d(mean[curves[[3]],],col=farbe[groups[i]],lwd=2)
#     spheres3d(mean[fix_drop,],radius=0.003, col=farbe[groups[i]])
#   }
# }
# 
# legend3d("topleft", levels(gp)[groups], pch = 16, col = farbe[groups], cex = 3)

# nearest neighbours based on procrustes distance
dist_mat <- as.matrix(dist(two.d.array(Proc$rotated[tokeep,,])))
rownames(dist_mat)
sort(dist_mat[which(rownames(dist_mat) == "UW105_811_LLM1"),])

# compare specimens
as.matrix(classifier$Name)
#specimens<-c(1,17); cbind(as.character(classifier$Name[c(specimens)]),farbe[c(1:length(specimens))])
specimens<-c(46,45); cbind(as.character(classifier$Name[c(specimens)]),farbe[c(1:length(specimens))])
show_landmarks = FALSE

clear3d("shapes");xmed<-median(Proc$rotated[,1,specimens[1]]); ymed<-median(Proc$rotated[,2,specimens[1]]); zmed<-median(Proc$rotated[,3,specimens[1]]);
plot3d(Proc$rotated[,,specimens[1]],type="n",xlab="", ylab="", zlab="",aspect =T,box=F,axes=F,xlim=c(xmed-0.4,xmed+0.4),ylim=c(ymed-0.4,ymed+0.4),zlim=c(zmed-0.4,zmed+0.4),col="red")

for(i in 1:length(specimens)){
  lapply(seq(curves), function(x){
    lines3d(Proc$rotated[curves_drop_orgn[[x]],,specimens[i]], col=farbe[i],lwd=2)
    if(show_landmarks){
      points3d(Proc$rotated[curves_drop_orgn[[x]],,specimens[i]],col=farbe[i], size = 7)
    }
    match_fixed <- curves_drop_orgn[[x]] %in% fix 
    if(any(match_fixed)){
      spheres3d(Proc$rotated[curves_drop_orgn[[x]][match_fixed],,specimens[i]],radius=0.003, col=farbe[i])
    }
  })
}
legend3d("topleft", cbind(as.character(classifier$Name[c(specimens)])), pch = 16, col = farbe, cex = 3)


#visualisation
a<-1;{
  clear3d("shapes");
  plot3d(Proc$rotated[curves_drop_orgn[[1]],,sub],aspect =F,col="black",type="l",xlab="",ylab="",zlab="",box=F,axes=F)
  lines3d(Proc$rotated[curves_drop_orgn[[2]],,sub],col="black")
  lines3d(Proc$rotated[c(curves[[1]],curves[[2]]),,sub],col="lightgray")
  lines3d(Proc$rotated[c(curves[[3]]),,sub],col="lightgray")
  spheres3d(SubCoordsArray[fix,,a],radius=0.05)
  texts3d(SubCoordsArray[unlist(curves_drop_orgn),,a],texts = unlist(curves_drop),font=2,cex=1)
  
  #texts3d(SubCoordsArray[unlist(curves_drop_orgn),,a],texts = unlist(curves_drop_orgn),font=2,cex=1)
}
