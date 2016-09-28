function bind_pe_procs()
  # filestream = open(ENV["PBS_NODEFILE"])
  home = ENV["HOME"]
  node_file_name = ENV["PE_HOSTFILE"]

  # parse the file - extract addresses and number of procs
  # on each
  # filestream = open("pe_file_ex.txt")
  filestream = open(node_file_name)
  seekstart(filestream)
  linearray = readlines(filestream)
  procs = map(linearray) do line
      line_parts = split(line," ")
      proc = {"name" => line_parts[1], "n" => line_parts[2]}
  end

  println("name of compute nodes and number of workers:")
  println(procs)

  # repeat for nodes with multiple procs
  # remove master from the node list
  master_node = ENV["HOSTNAME"]
  remove_master = 1
  machines = String[]
  for pp in procs
    # println(pp["name"])
    for i=1:int(pp["n"])
      if ( !contains(pp["name"],master_node)) | (remove_master==0)
        push!(machines,pp["name"])
      else
        remove_master=0
      end
    end
  end

  println("individual processes in machine file:")
  println(machines)


  println("adding machines to current system")
  addprocs(machines, dir= "/cm/shared/apps/Julia/0.3.0.pre/bin")
  println("done")
end

println("Started julia")

bind_pe_procs()


# here a function that runs your estimation:
# using MOpt, mig
include("examples/example-slices-mpi.jl")

println("done. quitting cluster.")

quit()