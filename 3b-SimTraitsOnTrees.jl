# get all the things you're using: 
using PhyloNetworks
using PhyloTraits
using DataFrames
using CSV
using StableRNGs

# read in your file with your shit and make a vector: 
dat = CSV.read("TreeTopologies_forSim.csv", DataFrame)
# get a STABLE RNG for this whole thing: 
rng = StableRNG(18480224)
# create your data frame for all of your node nums and names: 
NameNodeGuide = DataFrame()
ALL_TRAITS = DataFrame()
# set a rate vector (which will equal the number of simulations per tree, evenly divided by rate value.)
rates = repeat([0.1, 0.5, 1, 1.5, 2], outer=10)

# start a loop to go over each tree in dat: 
for d in 1:nrow(dat)
    # read your network: 
    net = readnewick(dat.TreeTopology[d])
# create your node number - name guide: 
     y = Vector{String}()
       
        for n in net.node
            name=n.name
            push!(y, name)
            end
# do the same thing for your node numbers: 
    x = Vector{Int64}()

        for n in net.node
            num = n.number
            push!(x, num)
            end
## combine them: 
    Names = DataFrame(NodeNum = x, NodeLabel = y)
## add your tree index for the node names and labels: 
    Names.Netname = repeat([dat.Netname[d]], outer= nrow(Names))
    # add this to your big df: 
    append!(NameNodeGuide, Names)
    # create an empty dataframe for your traits (this is re-made for each tree set): 
    Traits_big = DataFrame()
    

        for i in 1:50
            Params = ParamsBM( 0, rates[i])
            sim1 = rand(net, Params)

            #get your internal nodes: 
            int= DataFrame(trait=sim1[:internalnodes], nodeNum=sim1.M.internalnodenumbers)

            #get your tips: 
            tips = DataFrame(trait=sim1[:tips], nodeNum=sim1.M.tipnumbers)

            ## NOW COMBINE THESE
            TRAITS = [int;tips]
            TRAITS.Netname = repeat([dat.Netname[d]], outer= nrow(Names))
            TRAITS.SimNum = repeat([i], outer=nrow(TRAITS))
            TRAITS.rate = repeat([rates[i]], outer=nrow(TRAITS))
            append!(Traits_big, TRAITS)
        end
    # combine the 50 you just did with the rest of the nets:   
    append!(ALL_TRAITS, Traits_big)
 end

 # write your csv files: 
CSV.write("NodeLabelsNumbers_TreeSimsJulia.csv", NameNodeGuide)
CSV.write("TraitVals_TreeSims_Julia.csv", ALL_TRAITS)
