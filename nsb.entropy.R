################################################################################
# the functions below provide a simple R interface 
# to the C++ program "nsb-entropy" by Ilya Nemenman
# 
# written 2006-2008 by Jean Hausser 
# complains and insults should go to <jean@hausser.org>
#
# licence: GPL 3 (or later)
################################################################################

################################################################################
# usage:
#  source("nsb.entropy.R")
#  y = c(4, 2, 3, 0, 2, 4, 0, 0, 2) # bin counts
#  nsb.entropy(y)                   # 1.940646
################################################################################


##################### public function ##########################################

# input:   y     vector of counts (including zero counts)
#          CMD   path to nsb-entropy executable
# output:  nsb entropy estimate

# example usage:
#     y = c(4, 2, 3, 0, 2, 4, 0, 0, 2) # counts for each bin (including zeros!)
#     nsb.entropy(y)

nsb.entropy = function(y, CMD="./nsb-entropy")
{
  nsboutfile = nsbsave(file="nsbsamples.txt", y)
  system(paste(CMD, "-dnum -s1 -e1 -cY -v0 nsbsamples"))
  NSB.estimate = nsbload(file=nsboutfile)
  unlink(nsboutfile)
  unlink("nsbsamples.txt")
  return(NSB.estimate)
}
##################### private function #########################################

# saves a vector of counts (theta) such that it can be read by nsb-entropy
nsbsave = function(filename="samples.txt", theta) 
{
  datastr = character()
  for (i in 1:length(theta))
  {
    if ( theta[i] == 0 ) next
    for (j in 1:theta[i])
    {
      datastr = paste(datastr, i-1)
    }
  }
  #Removing leading space
  datastr = substr(datastr, 2, nchar(datastr))
  
  fileHandle = file(filename, "w")
  cat(
    paste(
      "# name: ASize",
      paste("# type: scalar", length(theta)),
      "# name: data",
      "# type: matrix",
      "# rows: 1",
      paste("# columns:", sum(theta)),
      datastr, 
      sep="\n"),
    file=fileHandle)
  close(fileHandle)

  # return the name nsb-entropy will store the results in
  nsboutfile = character()  
  splitname = strsplit(filename,'.',fixed=T)[[1]]
  return(paste(splitname[1], "_uni_num", length(theta),
                 "_mf1f0_1_entr.", splitname[2], sep=''))
}

# read text file output by nsb-entropy
nsbload = function(filename) 
{
  nsbout = readLines(filename)
  # Seek to the Snsb section
  for (i in 1:length(nsbout)) {
    if ( nsbout[i] == "# name: Snsb" ) break
  }
  #Snsb is at i+4
  return(as.double(nsbout[i+4]))
}
