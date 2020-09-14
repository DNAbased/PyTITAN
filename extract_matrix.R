# get arguments
args = commandArgs(TRUE)

# install relevant packages (as well as dependencies, if necessary)
if(!'TFBSTools' %in% installed.packages()){
    if(!'BiocManager' %in% installed.packages()) {
        install.packages('BiocManager', repos='https://cloud.r-project.org')
    }
    BiocManager::install('TFBSTools')
}
if(!'JASPAR2020' %in% installed.packages()) {
    if(!'remotes' %in% installed.packages()) {
        install.packages('remotes', repos='https://cloud.r-project.org')
    }
    remotes::install_github('da-bar/JASPAR2020')
}

# load packages
library(TFBSTools)
library(JASPAR2020)

# initialize jaspar datebase for subsequent use
db = 'data/matrixDB.sqlite'
initializeJASPARDB(db)

# only get human tfbs from the jaspar db
pwm_options = list('species' = 9606, 'all_versions' = TRUE, 'matrixtype' = toupper(args[1]))
pwm_list = getMatrixSet(JASPAR2020, pwm_options)

# write every tfbs matrix into a separate file for further processing in python
# names containing '::' do not work
for(list in pwm_list@listData) {
  write.table(data.frame(list@profileMatrix), paste0('data/', tolower(args[1]), '_format/', list@ID, '_', gsub('::', '.', list@name), '.pwm'), quote=FALSE, sep='\t', col.names=FALSE)
}
