# NetworksPhyloCompMethods
Here, you'll be able to find all code and data asosciated with the publication: NAME OF PUB. 
# What we did: 
We started by simulating a bunch of phylogenetic networks with different parameters with the R package SiPhyNetwork. For a detailed walkthrough of how to use and customize this program, see: https://github.com/jjustison/SiPhyNetwork, then we used those networks to simulate continuoulsy evolving traits under brownian motion in Julia package PhyloNetwworks. We extracted the major tree for each of our networks, and used the traits at the tip states to estimate ancestral characters using traditional tree-based PCMs. We used these data to understand how various features of trait evolutionary history and phylogenetic network structure impact PCM accuracy, and what we should do about it. 

If you want to follow along at every step, start by putting the two files "Net-trait-sims_ex1.txt.csv" and "Net-trait-sims_ex2.txt.csv" into a named directory. You'll need to be in those directories to follow scripts 4 and 6. 

# Simulating phylogenetic networks: 
1-SimNetsGetDat.R;
this script details how to simulate many networks at a time, collect important data from each network, collate that data across networks, and export both .net files, and relevant data files. 

# Simulating continuously evolving traits on networks: 
2-SimTraits_Net.jl;
This script details how to simulate continuous traits evolving on a phylogenetic network, extract relevant parameters and data for each simulation, and extract and output the major tree for each nework.

# Simulating continuously evolving traits on trees, and estimating ancestral character states: 
3-SimTraits_Tree.R;
This script details how to simulate continuoulsy evolving traits from a set of given topologies with a given set of parameters; and subsequently extract, and collate relevant information. Then, it uses that information to perform ancestral character estimations, comparing those to known parameters, and outputs a single file detailing true states, estimated states, error, (etc.) for each topology-trait combination. 

3b-SimTraitsOnTrees.jl; 
This does the same as above, but in phylonetworks on Julia instead of R.

*Note: This requires, as input, a dataframe in the form of: "TreeTopologies_ForSim.csv", in which one column contains the topology information, in newick text format, for all sampled topologies. 

# Estimating ancestral character states, rate parameters, and finding the best fit model of evolution: 
4-ACE_NETS.R;
This script details how to run ancestral character estimations in parallel for thousands of traits and topologies. It also compares those estimations to known parameters, and outputs a single file detailing true states, estimated states, error, (etc.) for each node and each topology-trait combination. 

*Note: This script requires input data in the form of: "Net-trait-sims_ex.csv", in which column names correspond to node and tip names, and values correspond to trait values, with additional information about simulation parameters, and one column containing the network topology in newick text format. There should be 1 .csv for network topology, containing all simulated traits for that single topology.  

# Combining network structure, trait history, and ACE data into a single file: 
5-Comb_NetDataAceData.R;
This script details how to combine data extracted from 1 and 4 into a single dataframe, containing all relevant information about network structure, trait history, and estimation error. 

# Calculating Node-Specific Variables: 
6-NodeSpecificVarCalc.R;
While step 5 collates data to analyze tree-wide and overall error across simulation scenarios, this script details how to use the output from steps 1 and 4 to collate a dataset for analysis of node-specific error in ancestral character estimation. 
*Note: This script uses the same files used in step 4. (i.e. "Net-trait-sims_ex.csv").

# Analyzing error distributions across datasets: 
7-CalcError_TreesandNets.R;
This script details how to take the output from 3 and 5, and use it to analyze erorr distributions in estimations of ancestral character states, evolutionary rate parameters, and model choice for ACEs using both tree-based and network-based simulations. 


# Important data: 
8-AllTopologies.csv;
A file containing all network topologies and corresponding major tree topologies used throughout this project. 

Net-trait-sims_ex1.txt.csv, Net-trait-sims_ex2.txtcsv, and TreeTopologies_ForSim.csv; 
These are sample datasets that you can download to follow along with each script. For the full datasets, see our Dryad page: (Link provided upon manuscript publication)




