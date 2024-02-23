
#PURPOSE:           Read NIIs and a predictor variable; perform mediation analysis with the voxel values as the mediator

#AUTHOR:            Seth D. Springer; modified original code of Camilo A Castelblanco & Lucas Weyrich
#VERSION HISTORY:   05/22/2022  v1: First working version



#Whole brain stats 
library(oro.nifti) #Necessary to read NIIs
library(openxlsx)
library(tibble)
library(tidyverse)
library(rstatix)
#library(profvis)
library(mediation)
library(bda)


setwd('D:')


nii_files = choose.files(default = "", caption = "Select NIIs for participants", multi = TRUE, filters = "*.nii")
data = read.csv(choose.files(default = "", caption = "Select CSV file with x and y values", multi = TRUE, filters = "*.csv"))



save_path = choose.dir(caption = 'Select the directory to save the statistical output')




#Effectively reading and loading one file
img <- readNIfTI(nii_files[1])
img_data = img@.Data
hdr = nifti_header(nii_files[1])
dimensions = dim(img@.Data) #obtain the dimensions nicely

n_participants = length(nii_files)

#Create a temp ID array for the ANOVA. Could use the actual GUID but this is more generalizable because we don't know GUID length 
tempID = 1:n_participants


#Preallocate 4D matrix
full_4D_matrix_c1 <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3], n_participants))


#Preallocate stats matrix

effect_coefficient_map <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3]))
p_value_map <- array(c(-1), dim = c(dimensions[1], dimensions[2], dimensions[3])) #can't be 0, otherwise can't threshold to see things with p-values
scaled_sigificance_map <- array(c(0), dim = c(dimensions[1], dimensions[2], dimensions[3]))


#test mediation

x = data$Age

y = data$group_correct_RT_OE




#Create a 4D matrix with brain data
for (i in 1:n_participants)
{
  #Assume that files are organized the same based on parID across conditions!! 
  
  #Load in pseudo-t for condition 1
  img_temp <- readNIfTI(nii_files[i])
  img_data_temp = img_temp@.Data
  
  full_4D_matrix_c1[ , , ,i] = img_data_temp
  
  rm(img_temp)
  
  
}


#Preallocate temp vectors
m = replicate(n_participants,0)



startTime <- Sys.time()

stats_counter = 1

pb <- winProgressBar(title = "Calculating voxel-wise statistics", min = 0, max = 1, initial = 0)

for (X_current in 1: dimensions[1]) 
{
  for (Y_current in 1: dimensions[2])
  {
    for (Z_current in 1: dimensions[3])
    {
      

      
      
      info <- sprintf("%d%% done", round(X_current/dimensions[1]))
      setWinProgressBar(pb, X_current/dimensions[1], label=info)
      
      
      #Extracting pseudo-t values from the current voxel
      m = full_4D_matrix_c1[X_current,Y_current,Z_current,]

      
      #To avoid model errors from zero padded voxels in the NIIs, just skip them
      #XX switch
      if (sum(m) == 0) { 
        
        #don't run stats for zero padded NII voxels
        
        
      } else  {
        
        
        
        
        #model.c         = lm(y ~  x)
        #model.c_summary = summary(model.c)
        
        
        model.a         = lm(m ~ x)
        model.a_summary = summary(model.a)
        model.a_X_pVal  = coef(model.a_summary)[2 , 4]
        
        model.Whole = lm(y ~ x + m)
        model.Whole_summary = summary(model.Whole)
        model.Whole_M_pVal  = coef(model.Whole_summary)[3 , 4]
        

        
        sobel_results = mediation.test(m,x,y)
        sobel_pVal = sobel_results$Sobel[2]
        
        
        
        #xx switch this back when done testing
        #if (model.a_X_pVal < .05 & model.Whole_M_pVal <.05) {
        #if (stats_counter>0) {
        if (sobel_pVal<.7) {
            
          
          #For testing
          #print(paste('X: ', X_current)) 
          #print(paste('Y: ', Y_current)) 
          #print(paste('Z: ', Z_current)) 
          
          
          
          
          
          stats_counter = stats_counter + 1 #counts the number of voxels that were bootstrapped
          
          
          results = mediate(model.a, model.Whole, treat='x', mediator='m', boot=TRUE, sims = 1000)
          #print(summary(results))
          
          

          p_value_map[X_current, Y_current, Z_current] =  results$d0.p
          
          
          #This scaled significance map is probably the most useful in neuroelf
          #If you only want p-values < 0.05, set the lower threshold to 1/0.05 = 20
          #If you only want p-values < 0.001, set the lower threshold to 1/0.001 = 1000
          #etc
          scaled_sigificance_map[X_current, Y_current, Z_current] = 1/results$d0.p
          
          
          
          effect_coefficient_map[X_current, Y_current, Z_current] = results$d0
          
          
        }
        
        

      }
    }
  }
}

#Close the loading bar
close(pb)




if (stats_counter == 1) {
  
print("Sorry but no mediation bootstrapping was ran since no voxels met the mediation requirements")

} else {

  
  
  #Write out the maps
  pixdim = pixdim(img)
  img.nifti = as.nifti(from = p_value_map, value = hdr,  verbose = FALSE) 
  img.nifti@reoriented = T
  img.nifti@datatype = 64
  img.nifti@bitpix = 64
  fname = file.path(path = save_path, 'WholeBrainMediation_pValue_wo_29_70_116_126')
  writeNIfTI(nim =img.nifti, fname, gzipped = FALSE) #gzipped = false to prevent the .nz in the file
  
  pixdim = pixdim(img)
  img.nifti = as.nifti(from = effect_coefficient_map, value = hdr,  verbose = FALSE) 
  img.nifti@reoriented = T
  img.nifti@datatype = 64
  img.nifti@bitpix = 64
  fname = file.path(path = save_path, 'WholeBrainMediation_effectCoeff_wo_29_70_116_126')
  writeNIfTI(nim =img.nifti, fname, gzipped = FALSE) #gzipped = false to prevent the .nz in the file  
  
  pixdim = pixdim(img)
  img.nifti = as.nifti(from = scaled_sigificance_map, value = hdr,  verbose = FALSE) 
  img.nifti@reoriented = T
  img.nifti@datatype = 64
  img.nifti@bitpix = 64
  fname = file.path(path = save_path, 'WholeBrainMediation_scaledSignificance_largerValueMeansMoreSig_wo_29_70_116_126')
  writeNIfTI(nim =img.nifti, fname, gzipped = FALSE) #gzipped = false to prevent the .nz in the file
  
}




endTime <- Sys.time()

total_time = endTime - startTime

print(total_time)



