
#PURPOSE:           Read NIIs and a predictor variable; perform regression and write out t and F stat maps

#AUTHOR:            Seth D. Springer; modified original code of Camilo A Castelblanco & Lucas Weyrich
#VERSION HISTORY:   05/22/2022  v1: First working version



#Whole brain stats 
library(oro.nifti) #Necessary to read NIIs
library(openxlsx)



all_participant_files = choose.files(default = "", caption = "Select NIIs for all participants", multi = TRUE, filters = "*.nii")
predictor_xlsx = choose.files(default = "", caption = "Select the excel with the predictor variable", multi = TRUE, filters = "*.xlsx")
covariate_xlsx = choose.files(default = "", caption = "Select the excel with the covariate variable", multi = TRUE, filters = "*.xlsx")


save_path = choose.dir(caption = 'Select the directory to save the statistical output')

data_predictor <- read.xlsx(predictor_xlsx, colNames = FALSE)
data_predictor = data_predictor$X1

data_covariate <- read.xlsx(covariate_xlsx, colNames = FALSE)
data_covariate = data_covariate$X1



#Effectively reading and loading one file
img <- readNIfTI(all_participant_files[1])
img_data = img@.Data
hdr = nifti_header(all_participant_files[1])
dimensions = dim(img@.Data) #obtain the dimensions nicely







n_participants = length(all_participant_files)


#Preallocate 4D matrix
full_4D_matrix <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3], n_participants))

#Preallocate stats matrix
t_stat_map <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3]))
F_stat_map <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3]))

#Create a 4D matrix with brain data
for (i in 1:n_participants)
{
  
  img_temp <- readNIfTI(all_participant_files[i])
  img_data_temp = img_temp@.Data
  
  full_4D_matrix[ , , ,i] = img_data_temp
  
  rm(img_temp)
  
}


#Preallocate temp vector
pseudo_Ts_temp = replicate(n_participants,0)

stats_counter = 1


pb <- winProgressBar(title = "Progress bar", min = 0, max = 1, initial = 0)

for (X_current in 1: dimensions[1]) 
{
  for (Y_current in 1: dimensions[2])
  {
    for (Z_current in 1: dimensions[3])
    {
      info <- sprintf("%d%% done", round(X_current/dimensions[1]))
      setWinProgressBar(pb, X_current/dimensions[1], label=info)
      
      pseudo_Ts_temp = full_4D_matrix[X_current,Y_current,Z_current,]
      
      #To avoid NANs and make more effecient!
      if (sum(pseudo_Ts_temp) == 0) { #I should be embarrased for forgetting the == !
         
        #dont do anything, leave the stats values at 0!
        
      } else  {
      model.brain_model = lm(pseudo_Ts_temp ~ data_predictor + data_covariate) #the number indicates the column!
      #model.brain_model = lm(empty_array[X_current, Y_current, Z_current, current_subject] ~ ages[current_subject])Syntax does not allow to have brackets inside
      t_value_temp = summary(model.brain_model)$coefficients["data_predictor","t value"]
      
      t_stat_map[X_current, Y_current, Z_current] = t_value_temp
      
      F_stat_map[X_current, Y_current, Z_current] = t_value_temp^2
      
      
      #only once, pull the degrees of freedom
      if (stats_counter ==  1) {
        df_1 = 1 #This is always one because we are interested in only the one predictor
        df_2 = as.numeric(summary(model.brain_model)$fstatistic[3])
        
        stats_counter = stats_counter+1
      }
      
      rm(model.brain_model) 
      }
    }
  }
}

close(pb)


t_stat_output_name = paste('lm_regression_with_covariate_t_values_', '_df_', df_2,sep = '')


#Write out the t stats brain
pixdim = pixdim(img)
img.nifti = as.nifti(from = t_stat_map, value = hdr,  verbose = FALSE) 
img.nifti@reoriented = T
fname = file.path(path = save_path, t_stat_output_name)
writeNIfTI(nim =img.nifti, fname, gzipped = FALSE) #gzipped = false to prevent the .nz in the file


F_stat_output_name = paste('lm_regression_with_covariate_F_values_', '_df_', df_1, '_',df_2,sep = '')


#Write out the t stats brain
pixdim = pixdim(img)
img.nifti = as.nifti(from = F_stat_map, value = hdr,  verbose = FALSE) 
img.nifti@reoriented = T
fname = file.path(path = save_path, F_stat_output_name)
writeNIfTI(nim =img.nifti, fname, gzipped = FALSE) #gzipped = false to prevent the .nz in the file

