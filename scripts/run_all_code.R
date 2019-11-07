#This is test code where each script is run and the time it takes to run is noted. There will 
#be more detailed testing (see Vignette or https://healtheconomicshackathon.github.io/hermes6/articles/Benchmarking.html)
#but the tests there are run multiple times to derive an average and takes longer to run.

# To run this it needs to be created as a package locally. 

# These lines 
if (!require("hermes6")) install.packages(".", repos = NULL, type="source")
if (!require("microbenchmark")) install.packages("microbenchmark")

# Restart your R session manually by going to Session and then Restart R but this code
# will do that.

.rs.restartR()

# Load package

library(hermes6)

# Creates list of functions where the name starts with markov. All top level functions should
# start with this but notes are needed for Vignette to ensure any no confusing results.

function_names <- ls("package:hermes6", pattern = "markov*")

# This function uses lapply which refers to above list and function that is called 'in line'
# func. Results are both printed and output to a list object.

output <- lapply(function_names, function(func) {
  func_call <- call(func)
  print(func_call)
  print(system.time(eval(func_call)))
})

