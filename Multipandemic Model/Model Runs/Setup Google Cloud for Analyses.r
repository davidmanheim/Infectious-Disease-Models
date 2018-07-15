#Setup for Google Cloud

Sys.setenv(GCE_AUTH_FILE="Authorization File.json")
Sys.setenv(GCE_DEFAULT_PROJECT_ID='Google Cloud Project ID')
Sys.setenv(GCE_DEFAULT_ZONE="us-east1-b")

library(testthat)
library(googleComputeEngineR)

#Initial Setup:
# vm <- gce_vm(template = "rstudio",
#              name = "my-rstudio",
#              username = "PickAUsername", password = "PickThePassword",
#              predefined_type = "f1-micro")

# Initial IP - Note the IP Address

# Load after that first time:
vm <- gce_vm(name = "my-rstudio")

#Stop and Start:
gce_vm_stop(vm)

gce_vm_start(vm)

#Initial setup of high-power machine: 

# bigvm <- gce_vm(template = "rstudio",
#              name = "my-big-rstudio",
#              username = "PickAUsername", password = "PickThePassword",
#              predefined_type = "n1-highcpu-2")
# IP - Note the IP Address

#Stop It!
gce_vm_stop(bigvm)

#Restart:
bigvm <- gce_vm(name = "my-big-rstudio")
# Note the IP Address