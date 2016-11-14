
universe = vanilla
getenv   = true

machine_count  = 1
request_cpus   = 1
request_memory = 1GB

executable = main_NoReplace_$(Process)

log    = NoReplace_$(Process).log 
output = NoReplace_$(Process).out
error  = NoReplace_$(Process).err

queue 10
