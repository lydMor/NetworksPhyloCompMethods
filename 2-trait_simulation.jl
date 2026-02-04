#written by Sungsik Kong 2024

using PhyloNetworks
using Distributions
using DataFrames
using Random
using Dates
using CSV

function f(;
            foldername=["BLG", "BLN", "BLN75", "LG75", "LN50", "BLG75", "BLN50", "LG", "LN", "LN75"],
            numreps=100::Integer, 
            ancmean=[0]::Array, 
            variances=[0.1, 0.5, 1, 1.5, 2]::Array, 
            proptransgressive=0.1::Float64,
            #transgmin=(-10.0)::Float64,
            #transgmax=10.0::Float64,
            wdir="/Users/khaosan/Dropbox/sungsik.kong-UWisc/2023-Lydia, Emma, Kevin-npcm/sim_traits/",
            csv=false::Bool
            )

    #move to the working directory that contains folders named in the variable foldername
    for folder in foldername
        #move into the folder
        cd(wdir)
        cd(folder)
        
        #create a logfile and record basic information
        logfile = open("../res/$(folder).log", "w+")

        log="Timestamp: $(Dates.format(now(), "HH:MM U d Y"))\n"; pwf(log,logfile)
        log="[setting wdir] Working directory: $(pwd())\n"; pwf(log,logfile)
        log="[setting numreps] number of replicates: $numreps\n"; pwf(log,logfile)
        log="[setting variances] variances: $variances\n"; pwf(log,logfile)
        log="[setting proptransgressive] probability of transgressive evolution: $proptransgressive\n"; pwf(log,logfile)
        log="[setting csv] will we store the result in a csv file?: $csv\n"; pwf(log,logfile)
        log="We recommend to set csv=true as not a portion of long newick won't be visible in STDOUT.\n"; pwf(log,logfile)
        #list the network files in the folder (plus many other unwanted residual files)
        networks=readdir()
        
        log="List of files in directory: $(networks)\n"; wf(log,logfile)
        log="fyi: There are $(length(networks)) files in this working directory.\n"; pwf(log,logfile)
        delim="\n"; pwf(delim,logfile)

        for network in networks  
            log="Timestamp: $(Dates.format(now(), "HH:MM U d Y"))\n"; pwf(log,logfile)
            log="Continuous trait simulation using the network in $(network)\n"; pwf(log,logfile)
            
            #read in the network
            log="Reading network in $network..."; 
            #net=readTopology(network)
            
            try readTopology(network)
            catch err 
                msg="Error while reading the file $(network). May be it does not contain a newick or corrupted. Moving onto the next file.\n"; pwf(msg,logfile)
                log*="failed\n"; pwf(log,logfile)
                delim="\n"; pwf(delim,logfile)
                continue
            end
            net=readTopology(network)
            log*="success\n"; pwf(log,logfile)
            
            #PhyloNetworks.printEverything(net)

            #by default the leaves and hybrid nodes are already labelled (lx for leaves, Hx for hybridnodes, where x=integer)
            #we assign the internal tree node names as i(abs(node.number))
            #we do abs to turn the negative node number to positive to renove potential complexity in downstream analysis
            function assignintnodelabels(net::HybridNetwork)
                for eachnode in net.node
                    if isempty(eachnode.name)
                        !(eachnode.leaf) || error("Node $(eachnode.number) is a leaf but unlabelled.")
                        !(eachnode.hybrid) || error("Node $(eachnode.number) is a hybrid node but unlabelled.")
                        eachnode.name="i$(abs(eachnode.number))"
                    else
                        continue
                    end
                end
            end
            #PhyloNetworks.printEverything(net)

            #assigning internal node names in the network
            log="Assigning internal tree node names in $network..."
            try assignintnodelabels(net)
                log*="success\n"; println(log); pwf(log,logfile)
            catch err 
                log*="failed\n"; println(log); pwf(log,logfile)
                println(err)
            end

            ###major tree
            t=PhyloNetworks.displayedTrees(net,0.0) #find num dts
            log="There are $(length(t)) displayed trees in the network $network.\n"; pwf(log,logfile)
            majT=PhyloNetworks.majorTree(net,nofuse=true) #nofuse to compute the gamma
            #get the prob of the major tree
            MajTreeprop=1
            for edge in majT.edge
                MajTreeprop*=edge.gamma
            end
            log="Major tree with gamma=$MajTreeprop is selected and recorded.\n"; pwf(log,logfile)
            majT=PhyloNetworks.majorTree(net,nofuse=false) #fuse to get the topology

            ###create dataFrame
            log="Creating an empty dataframe for simulation using $network..."
            df=DataFrame(filename=String[],
                         NetTopology=String[],
                         NumDisplayedTrees=Integer[],
                         MajTreeTopology=String[],
                         MajTreeProp=Float64[],
                         seed=Integer[],
                         ancestralmean=Float64[],
                         variance=Float64[],
                         replicate=Integer[],
                         numHybrids=Integer[],
                         numTransgressive=Integer[]
                        )
            
                        for eachnode in net.node
                            df[!,eachnode.name].=0.0
                            if eachnode.hybrid
                                df[!,"br$(eachnode.name)"].=0.0
                            else
                                continue
                            end
                        end
            log*="success\n"; pwf(log,logfile)
#=
In Jupyter notebook, show:
    show empty df
    explain iX, HX, brHX, etc.
=#

            #the information so far that needs to tbe go into the spreadsheet
            filename=network
            NetTopology=PhyloNetworks.writeTopology(net)
            NumDisplayedTrees=length(t)
            MajTreeTopology=PhyloNetworks.writeTopology(majT)
            MajTreeprop=MajTreeprop
            numHybrid=net.numHybrids

            #println(filename)
            #println(NetTopology)
            #println(NumDisplayedTrees)
            #println(MajTreeTopology)
            #println(MajTreeprop)

            #simulation 
            for mean in ancmean
                for variance in variances
                    transgmin=(-1)*variance*20
                    transgmax=variance*20
                    log="[setting transgmin] minimum value that can get by transgressive evol: $transgmin\n"; pwf(log,logfile)
                    log="[setting transgmax] maximum value that can get by transgressive evol: $transgmax\n"; pwf(log,logfile)
                    log="Continuous trait simulation using the network in $(network)\n"; pwf(log,logfile)
                    log="Simulating $numreps replicates of continuous trait data using BM with mean=$mean and variance=$variance.\n"; pwf(log,logfile)
                    for replicate in 1:numreps
                        log="Simulation replicate $replicate/$numreps initiated.\n"; pwf(log,logfile)
                        #set seed
                        seed=Int(rand(1:10e5))
                        Random.seed!(seed)
                        log="Setting seed to $seed.\n"; pwf(log,logfile)
                        #log="Seed: $seed"; println(log)
                        #set our dictionary of the node number and trait values
                        d = Dict{String,Float64}()

                        #transgressive evolution
                        log="Simulating transgressive evolution with prop=$(proptransgressive), minval=$(transgmin), and maxval=$(transgmax)..."
                        transgressive=zeros(numHybrid)
                        numtransgressive=0

                        for apotentialtransgressiveevent in 1:numHybrid
                            roll=rand(Uniform(0, 1))
                            if roll<=proptransgressive
                                numtransgressive+=1
                                shiftval=variance
                                while abs(shiftval) <= variance
                                    shiftval=rand(Uniform(transgmin, transgmax))
                                end
                                transgressive[apotentialtransgressiveevent]=shiftval
                            else
                                transgressive[apotentialtransgressiveevent]=0
                            end
                        end
                        log*="success\n"; pwf(log,logfile)

                        #for elem in transgressive
                        #    if !iszero(elem) numtransgressive += 1
                        #    else continue
                        #    end
                        #end
                        log="Total $numtransgressive/$(net.numHybrids) transgressive evolution occurred.\n"; pwf(log,logfile)

                        shift = PhyloNetworks.shiftHybrid(transgressive, net)
                        shiftedge=PhyloNetworks.getShiftEdgeNumber(shift)
                        shiftval=PhyloNetworks.getShiftValue(shift)
                        brnametrasgressivepairs=Tuple[]
                        #println(shiftedge)
                        length(shiftedge)==length(shiftval) || error("Something went wrong while simulating transgressive evolution.")
                        for i in 1:length(shiftval)
                            parent=PhyloNetworks.getParent(net.edge[shiftedge[i]])
                            push!(brnametrasgressivepairs,("br$(parent.name)",shiftval[i]))
                        end
                        
                        for tuple in brnametrasgressivepairs
                            push!(d,tuple[1]=>tuple[2])
                        end
                        #println(shiftval)
                        #println(length(shiftval))
                        #println(brnametrasgressivepairs)
                        #println(transgressive)
                        #println(numtransgressive)
                        #println("$(numtransgressive) (max: $numHybrid) transgressive evolution events occurred for this trait in this scenario.")
                        #display(shift)

                        #trait simulation with shifts
                        if !iszero(numtransgressive)
                            log="Continuous trait simulation with any character shift at a hybrid node..."
                        else
                            log="Continuous trait simulation without character shift..."
                        end

                        params = PhyloNetworks.ParamsBM(mean,variance,shift)
                        sim = PhyloNetworks.simulate(net,params)
                        log*="success\n"; pwf(log,logfile)
                        #PhyloNetworks.printEverything(net)
                        #println(sim.M.nodeNumbersTopOrder)
                        #println(tipLabels(sim))
                        #println(sim.M.tipNumbers)
                        #println(sim[:Tips])
                        #println(sim.M.internalNodeNumbers)
                        #println(sim[:InternalNodes])
                        numinternalnodes=net.numNodes-net.numTaxa
                        #println(numinternalnodes)
                        
                        #create a disctionary and store it somewhere 
                        nodenumtraitpairs=Tuple[]
                        
                        length(sim.M.tipNumbers)==(net.numTaxa) || error("Number of simulated trait at tips does not match with the number of tips.")
                        length(sim[:Tips])==(net.numTaxa) || error("Number of simulated trait at tips does not match with the number of tips.")
                        length(sim.M.internalNodeNumbers)==(numinternalnodes) || error("Number of simulated trait at tips does not match with the number of nonleaf nodes.")
                        length(sim[:InternalNodes])==(numinternalnodes) || error("Number of simulated trait at tips does not match with the number of nonleaf nodes.")
                        
                        #create a tuple for each (taxon name, trait value) pair and collect at nodenumtraitpairs
                        for i in 1:net.numTaxa
                            for nod in net.node
                                if nod.number==sim.M.tipNumbers[i]
                                    push!(nodenumtraitpairs,(nod.name,sim[:Tips][i]))
                                else
                                    continue
                                end
                            end
                        end
                        #create a tuple for each (internal node name, trait value) pair and collect at nodenumtraitpairs
                        for i in 1:numinternalnodes
                            for nod in net.node
                                if nod.number==sim.M.internalNodeNumbers[i]
                                    push!(nodenumtraitpairs,(nod.name,sim[:InternalNodes][i]))
                                else
                                    continue
                                end
                            end
                        end

                        #push to the dictionary
                        for tuple in nodenumtraitpairs
                            push!(d,tuple[1]=>tuple[2])
                        end
                        #println(d)
                        
                        #record results
                        data=[filename,
                                NetTopology,
                                NumDisplayedTrees,
                                MajTreeTopology,
                                MajTreeprop,
                                seed,
                                mean,
                                variance,
                                replicate,
                                numHybrid,
                                numtransgressive#11
                                ]

                        dfhdr=names(df)
                        for i in 12:length(dfhdr)
                            try push!(data,d[dfhdr[i]]) #push simulated trait to df for the cooresponding node
                            catch err push!(data,0.0) 
                            #if no node number found in dic, this tells us the branch following the hybrid node
                            #where transgressive evolution did not occur, thus we put value of shift as zero
                            end
                        end
                        
                        #println(data)
                        push!(df,data)
                    end
                    log="Simulation...complete\n"; pwf(log,logfile)
                    delim="\n"; pwf(delim,logfile)
                end
            end

            if !(csv)
                res="$df\n"; pwf(res,logfile)
            else
                res="$df\n"; pwf(res,logfile)
                CSV.write("$network.csv", df)
                println("CSV file written in the current directory with the filename $network.csv\n")
            end
            delim="\n"; pwf(delim,logfile)
            #println("----------------------------------------\n")
            
        end
        close(logfile)    
    end
end


#print, write, and flush
function pwf(log::String,logfile::IOStream)
    print(log); write(logfile, log); flush(logfile)
end

#write and flush but no print STDOUT
function wf(log::String,logfile::IOStream)
    write(logfile, log); flush(logfile)
end