
#PURPOSE:           Read NIIs and a predictor variable; perform regression and write out t and F stat maps

#AUTHOR:            Seth D. Springer; modified original code of Camilo A Castelblanco & Lucas Weyrich
#VERSION HISTORY:   05/22/2022  v1: First working version



#Whole brain stats 
library(oro.nifti) #Necessary to read NIIs
library(openxlsx)
library(tibble)
library(tidyverse)
library(rstatix)
#library(profvis)


condition1_files = choose.files(default = "", caption = "Select NIIs for condition 1 participants", multi = TRUE, filters = "*.nii")
condition2_files = choose.files(default = "", caption = "Select NIIs for condition 2 participants", multi = TRUE, filters = "*.nii")
condition3_files = choose.files(default = "", caption = "Select NIIs for condition 3 participants", multi = TRUE, filters = "*.nii")

save_path = choose.dir(caption = 'Select the directory to save the statistical output')


#Effectively reading and loading one file
img <- readNIfTI(condition1_files[1])
img_data = img@.Data
hdr = nifti_header(condition1_files[1])
dimensions = dim(img@.Data) #obtain the dimensions nicely

n_participants = length(condition1_files)

#Create a temp ID array for the ANOVA. Could use the actual GUID but this is more generalizable because we don't know GUID length 
tempID = 1:n_participants


#Preallocate 4D matrix
full_4D_matrix_c1 <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3], n_participants))
full_4D_matrix_c2 <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3], n_participants))
full_4D_matrix_c3 <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3], n_participants))

#Preallocate stats matrix

F_stat_map <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3]))

#Create a 4D matrix with brain data
for (i in 1:n_participants)
{
  #Assume that files are organized the same based on parID across conditions!! 
  
  #Load in pseudo-t for condition 1
  img_temp <- readNIfTI(condition1_files[i])
  img_data_temp = img_temp@.Data
  
  full_4D_matrix_c1[ , , ,i] = img_data_temp
  
  rm(img_temp)
  
  #Load in pseudo-t for condition 2
  img_temp2 <- readNIfTI(condition2_files[i])
  img_data_temp2 = img_temp2@.Data
  
  full_4D_matrix_c2[ , , ,i] = img_data_temp2
  
  rm(img_temp2)
  
  #Load in pseudo-t for condition 3
  img_temp3 <- readNIfTI(condition3_files[i])
  img_data_temp3 = img_temp3@.Data
  
  full_4D_matrix_c3[ , , ,i] = img_data_temp3
  
  rm(img_temp3)
  
}


#Preallocate temp vectors
pseudo_Ts_temp1 = replicate(n_participants,0)
pseudo_Ts_temp2 = replicate(n_participants,0)
pseudo_Ts_temp3 = replicate(n_participants,0)

stats_counter = 1

pb <- winProgressBar(title = "Calculating voxel-wise statistics", min = 0, max = 1, initial = 0)

for (X_current in 1: dimensions[1]) 
{
  for (Y_current in 1: dimensions[2])
  {
    for (Z_current in 1: dimensions[3])
    {
      
      #(paste('X: ', X_current)) 
      #print(paste('Y: ', Y_current)) 
      #print(paste('Z: ', Z_current)) 
      
      
      info <- sprintf("%d%% done", round(X_current/dimensions[1]))
      setWinProgressBar(pb, X_current/dimensions[1], label=info)
      
      pseudo_Ts_temp1 = full_4D_matrix_c1[X_current,Y_current,Z_current,]
      pseudo_Ts_temp2 = full_4D_matrix_c2[X_current,Y_current,Z_current,]
      pseudo_Ts_temp3 = full_4D_matrix_c3[X_current,Y_current,Z_current,]
      
      #To avoid NANs and make more effecient!
      if (sum(pseudo_Ts_temp1) == 0 || sum(pseudo_Ts_temp2) == 0 || sum(pseudo_Ts_temp3) == 0) { #I should be embarrased for forgetting the == !
        
        #dont do anything, leave the stats values at 0!
        
      } else  {
        
        df = data.frame(tempID,pseudo_Ts_temp1,pseudo_Ts_temp2,pseudo_Ts_temp3)
        #df = as_tibble(df)
        
        df <- df %>%
          gather(key = "condition", value = "pseudoT", pseudo_Ts_temp1, pseudo_Ts_temp2, pseudo_Ts_temp3) %>%
          convert_as_factor(tempID, condition)
        
        #res.aov <- anova_test(data = df, dv = pseudoT, wid = tempID, within = condition)
        
        model <- aov(pseudoT~factor(condition)+Error(factor(tempID)), data = df)
        
        #view model summary
        model_summary = summary(model)
        model_summary = model_summary[[2]]
        
        F_value = model_summary[[1]][1,4]
        
        
        rm(df)
        
        
        F_stat_map[X_current, Y_current, Z_current] = F_value
        
        #only once, pull the degrees of freedom
        if (stats_counter ==  1) {
          df_1 = model_summary[[1]][1,1]
          df_2 = model_summary[[1]][2,1]
          
          stats_counter = stats_counter+1
        }
        rm(model)
      }
    }
  }
}

close(pb)


F_stat_output_name = paste('ANOVA_RepeatedMeasures_F_values_', '_df_', df_1, '_',df_2,sep = '')

#Write out the F stats brain
pixdim = pixdim(img)
img.nifti = as.nifti(from = F_stat_map, value = hdr,  verbose = FALSE) 
img.nifti@reoriented = T
fname = file.path(path = save_path, F_stat_output_name)
writeNIfTI(nim =img.nifti, fname, gzipped = FALSE) #gzipped = false to prevent the .nz in the file

