library('NuPoP')
nupop_input_dir <- '../nupop_input_temp/'

# list all files in the input directory
files <- list.files(nupop_input_dir, pattern = "*.seq")

# for loop to run NuPoP on each file
for (file in files) {
    file_ = paste0(nupop_input_dir,file)
  # run predNuPoP on each file
    predNuPoP(file_,species=1,model=4)
}

# print name of each file
for (file in files) {
    file_ = paste0(nupop_input_dir,file)
    print(file_)
}