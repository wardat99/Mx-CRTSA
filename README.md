# Mx-CRTSA
The source code for our paper "Discovering Community Structure in Multiplex Networks via a Co-Regularized Robust Tensor-Based Spectral Approach".

Available [here](https://www.mdpi.com/2076-3417/13/4/2514).

# Required toolboxes:
The MATLAB code for MV-RTSC needs some toolboxes:

[1] Tensor decomposition, Xie, Y., Tao, D., Zhang, W., Liu, Y., Zhang, L., & Qu, Y. (2018). 
"On unifying multi-view self-representations for clustering by tensor multi-rank minimization."

[2] Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

# Data set:
All datasets used in our paper are available at Google Drive. Each dataset is a mat file containing 2 variables including adgacancy tensor and ground truth of the network. https://drive.google.com/drive/folders/1T8AITuQZCbCB52PELprIVb75d4-6hi20


# Citation:

If you find this code helpful, please cite the paper as follows:

```

@Article{app13042514,
AUTHOR = {Al-sharoa, Esraa and Al-wardat, Mohammad and Al-khassaweneh, Mahmood and Bataineh, Ali Al},
TITLE = {Discovering Community Structure in Multiplex Networks via a Co-Regularized Robust Tensor-Based Spectral Approach},
JOURNAL = {Applied Sciences},
VOLUME = {13},
YEAR = {2023},
NUMBER = {4},
ARTICLE-NUMBER = {2514},
URL = {https://www.mdpi.com/2076-3417/13/4/2514},
ISSN = {2076-3417},
ABSTRACT = {Complex networks arise in various fields, such as biology, sociology and communication, to model interactions among entities. Entities in many real-world systems exhibit different types of interactions, which requires modeling these type of systems properly. Multiplex networks are used to model these systems, as they can reflect the nodes&rsquo; pair-wise interactions as multiple distinct types of links across layers. Community detection is a widely studied application in network analysis as it provides insights into the structure and organization of the network. Even though multiple algorithms have been developed in the community detection field, many of them have a limited performance in the presence of noise. In this article, we develop a novel algorithm that combines tensor low-rank representation, spectral clustering and distance regularization to improve the accuracy in discovering communities in multiplex networks. The low-rank representation leads to reducing the noise and errors existing in the network and the optimization of an accurate consensus set of eigenvectors that reveals the communities in the network. Moreover, the proposed approach balances the agreement between the eigenvectors of each layer, i.e., individual subspaces, and the consensus set of eigenvectors, i.e., common subspaces, by minimizing the projection distance between them. The common and individual subspaces are computed efficiently through Tucker decomposition and modified spectral clustering, respectively. Finally, multiple experiments are conducted on real and simulated networks to evaluate the proposed approach and compare it to state-of-the-art algorithms. The proposed approach shows its robustness and efficiency in discovering the communities in multiplex networks.},
DOI = {10.3390/app13042514}
}
}
```
