# Code for analyzing ecological multilayer networks

This repository contains a (hopefully growing) set of algorithms to analyze multilayer ecological networks.

## NEE 2017

This is a collection of R code and of Matlab functions used to analyze the networks in the publication: *[The Multilayer Nature of Ecological Networks](https://www.nature.com/articles/s41559-017-0101)*. Predominantly this is for the analysis of community detection (modularity). But also there is code for robustness analysis and reducibility. If you want to replicate the results in the paper you should use the code with the data set in the folder **Data**

Please read the original paper **AND IT'S SUPPLEMENTARY INFORMATION** for details on what these functions do. Also, each code has a detailed description but if you need more details or get stuck please open an issue. This will also help others.

The code which creates the modularity matrix for bipartite networks was modified from http://netwiki.amath.unc.edu/GenLouvain/GenLouvain with help of Lucas Jeub and Mason A. Porter.

The GenLouvain code (version 2.0, 2014) was downloaded from: http://netwiki.amath.unc.edu/GenLouvain/GenLouvain on July 2016 and was not modified. Please cite this code when you use it (or parts of it). Citation:

Inderjit S. Jutla, Lucas G. S. Jeub, and Peter J. Mucha, "A generalized Louvain method for community detection implemented in MATLAB," http://netwiki.amath.unc.edu/GenLouvain (2011-2014).

The code for **reducibility** was adapted from Manlio De Domenico's excellent muxViz software, and based on the paper:
De Domenico, M., Nicosia, V., Arenas, A. & Latora, V. (2015). Structural reducibility of multilayer networks. Nat. Commun., 6:6864, doi:10.1038/ncomms7864


Please cite the NEE paper as well as the other relevant paper on which it is based when you use the code. Citation:

 *[The Multilayer Nature of Ecological Networks](https://www.nature.com/articles/s41559-017-0101) (2017) Pilosof S, Porter MA, Pascual M, Kefi S. Nature Ecology & Evolution 1:0101.* doi:10.1038/s41559-017-0101

For any questions/details pleae don't hesitate to contact me at: shainova@gmail.com
