# Codes for paper *Fibration symmetry uncovers minimal regulatory networks for logical computation in bacteria*

The flow is the following:
1. Get fibers for the network using the script **fiber.R** from the directory **Codes_for_fibers**
2. Use **Clean_paper.R** to obtain the *Minimal Network* (the core of the original entire network).
3. Use **circuits_luis.R**, **cycles.R** to find the circuits and the cycles in the network, correspondigly.
4. The **randomized.R** script is for obtaining a randomized version of the original network that still preserves the degree distribution (i.e. a re-wiring) to know what would be the core structures for randomly generated networks and thus compare the likelihood of the actual observed structures in the real networks arising by chance

## 1. Codes_for_fibers.R

This is a modification of an older version of the codes that were made available on [this](https://github.com/makselab/fibrationSymmetries) repository. Most of the changes made were so that it would be easier to work with the rest of the codes for the paper. (Its important to note that all the C++ scripts need to be on the same directory as **fiber.R**)
The script **fiber.R** runs all of the scripts in C++ to obtain the coloring of the nodes, or fibers, of the network. This script calls the functions on the **functions.R** script, firstly for preparing the network files to be run on the C++ codes and later to retrieve the obtained colors. After the colors are obtained the script calls **classifier.R** to classify the building blocks. 

## 2. Clear_paper.R

This script requires both the network file (in an edge -> edge format) and a file with the list of nodes and their colors or fibers. The script first collapses the network to its base, i.e. collapses all the fibers to a single representative node per fiber. This reduced network still preserves the same information flow dynamics. The base of the network is further reduced by applying the kcore decomposition to obtain the k<sub>out</sub>=0 core of the network. 

## 3. circuits_luis.R and cycles.R

Both of this scripts take the network file and returns a list of all the circuits/cycles found, it also classify the circuits according to the type of edges (different types of edges produce different logical circuits).

## 4. randomized.R 

Again, takes the network files, outputs an ensemble of random networks (a rewiring of the original network, i.e. the null model) that still preserve the degree distribution and then analyzes the structure at the core of these randomized networks to give a z-score as a notion of how close to the expected structure for a random network the observed structure is. 
