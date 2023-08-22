# Codes for paper *Fibration symmetry uncovers minimal regulatory networks for logical computation in bacteria*
---
The flow is the following:
1. Get fibers for the network using the script **fiber.R** from the directory **Codes_for_Fibers**
2. Use **Clean_paper.R** to obtain the *Minimal Network* (the core of the original entire network).
3. Use **circuits_luis.R**, **cycles.R** to find the circuits and the cycles in the network, correspondigly.
4. The **randomized.R** script is for obtaining a randomized version of the original network that still preserves the degree distribution (i.e. a re-wiring) to know what would be the core structures for randomly generated networks and thus compare the likelihood of the actual observed structures in the real networks.


[this](https://github.com/makselab/fibrationSymmetries)

