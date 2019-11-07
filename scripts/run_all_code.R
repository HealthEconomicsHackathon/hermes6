# This is test code where each script is run and the time it takes to run is noted. There will 
# be more detailed testing (see Vignette or https://healtheconomicshackathon.github.io/hermes6/articles/Benchmarking.html)
# but the tests there are run multiple times to derive an average and takes longer to run.

# To run this it needs to be created as a package locally. Note that if the entire script is 
# selected to run it will stop at the restartR() line of code.

# Install packages if not already installed.

if (!require("hermes6")) install.packages(".", repos = NULL, type = "source")
if (!require("tidyverse")) install.packages("tidyverse")

# Restart your R session manually by going to Session and then Restart R but this code
# will do that.

.rs.restartR()

# Load packages

library(hermes6)
library(tidyverse)

# Creates list of functions where the name starts with markov. All top level functions should
# start with this. 

function_names <- ls("package:hermes6", pattern = "markov*")

# This function uses lapply which refers to above list and function that is called 'in line'
# func. Results are both printed and output to a list object and only elapsed time is required.

output <- lapply(function_names, function(func) {
  func_call <- call(func)
  print(func_call)
  out <- system.time(eval(func_call))
  print(out)
  out["elapsed"]
})

# Names of the functions are not with the output above. This line attaches the function
# names to the list items.

names(output) <- function_names


# Export will be in wide form but by moving to a df using enframe() this goes to long form.

output_df <- output %>% 
  unlist(recursive = FALSE) %>% 
  enframe() 

# Export to .csv

write.csv(output_df, file = "RunningTimeModels.csv",row.names=TRUE)
