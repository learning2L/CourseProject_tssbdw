# Proposal

## Introduction



Clustering

DP 1. stick breaking process 2. used to cluster

HDP & nDP: 

* Important extension of DP!
* 也可以用于聚类
* 各自的drawback，而且都没办法给出 tree-structured clusters

Adams proposed tree-structured stick process (TSSB) prior：

* two-dimension extension of DP; tree-structured clusters
* drawback: slice sampler is slow; hard to handle

We consider finite truncation of TSSB, referred to as TSSB-DW (depth and width)

prove that the difference between TSSB and TSSB-DW has upper bound. $\checkmark$

New sampler method（block sampler） + HMC to sample parameters (HMC 的优点) <font color='blue'> variational inference ? </font>

Real data: 

* Clustering of images (CIFAR-10) --- 了解数据的收集过程
* 提一下 生物数据（更高维）——应用广泛
  * 维度更高（也可以筛选降维）
  * 不止是 gene expression data (sc-RNA)，也可以是 进化树的数据等等，都可以建立 tree-structured clusters



---

Dirichlet process (DP) is one of the most widely-used nonparametric Bayesian prior for solving clustering problems, since it is naturally  discrete. Due to this property, DP has its stick-breaking representation.

> There are two extensions of DP. Hierarchical Dirichlet process (HDP) was proposed by Teh et al. (2006). The distributions share the same atoms, but the weights are different from each other. Nested Dirchlet process (nDP), presented by Rodriguez et al. (2008), 

Alternatively, Adams et al. (2010) extended the stick-breaking process of DP to a two-dimension case, which is referred to as tree-structured stick-breaking process (TSSB). TSSB has infinite width and depth, and each datum can locate in any internal node. The authors proposed a slice sampler to carry out the posterior inference, and applied the model to image-hierarchical-clustering problem (CIFAR-100). Although TSSB can solve the tree-structured clustering problem, because of the infinity of depth and width, the slice sampling process might be slightly slow and difficult to implement.

Truncation of stick-breaking process (Ishwaran and James, 2001) provided an excellent direction to avoid the infinity property. Based on the truncation version, Ishwaran and James presented the block sampler, which can be further combined with parallel computing.

Motivated by (Ishwaran and James, 2001), we want to consider a truncation version of the depth and width of TSSB, which is referred to as TSSB-DW (truncation of **TSSB** in **D**epth and **W**idth), and design a more efficient sampling approach. 

We use the equivalent form of the model, where xxx are independent. So we can implement parallel computation to save time. 

We also use HMC to sample parameters, since HMC uses the gradient information of the posterior distribution, which is more efficient.

We use the CIFAR-10 dataset to test our model. Further, this model is also suitable for gene expression data whose dimension can be larger.





---

## EDA



<font color='red'>新的想法：用师兄paper中的数据，用TSSB跑一下试试！！！</font>



## Model



### TSSB

$G=\sum_{\epsilon} \pi_{\epsilon} \delta_{\theta_{\epsilon}},\quad \pi_{\epsilon}\sim \text{TSSB}(\alpha,\gamma)$ 

where $\pi_{\epsilon}=\nu_{\epsilon} \varphi_{\epsilon} \prod_{\boldsymbol{\epsilon}^{\prime} \prec \boldsymbol{\epsilon}} \varphi_{\boldsymbol{\epsilon}^{\prime}}\left(1-\nu_{\boldsymbol{\epsilon}^{\prime}}\right), \quad \varphi_{\boldsymbol{\epsilon} \epsilon_{i}}=\psi_{\boldsymbol{\epsilon} \epsilon_{i}} \prod_{j=1}^{\epsilon_{i}-1}\left(1-\psi_{\boldsymbol{\epsilon} j}\right), \quad \pi_{\emptyset}=\nu_{\emptyset}$

and $\theta_{\epsilon \epsilon_i}|\theta_{\epsilon} \sim p(\cdot|\theta_{\epsilon}), \quad \nu_{\boldsymbol{\epsilon}} \sim \operatorname{Be}(1, \alpha(|\boldsymbol{\epsilon}|)),\quad \psi_{\boldsymbol{\epsilon}} \sim \operatorname{Be}(1, \gamma)$  

### TSSB-DW

$G=\sum_{\epsilon} \pi_{\epsilon} \delta_{\theta_{\epsilon}},\quad \pi_{\epsilon}\sim \text{TSSB-DW}(\alpha,\gamma)$ 

where $\nu_{\boldsymbol{\epsilon}} \sim \operatorname{Be}(1, \alpha(|\boldsymbol{\epsilon}|)),\ |\epsilon|=0,1,...,D-1,\quad \nu_{\boldsymbol{\epsilon}}=1$ when $|\epsilon|=D$ 

and given $\epsilon$, $\psi_{\boldsymbol{\epsilon}\epsilon_i} \sim \operatorname{Be}(1, \gamma),\ \epsilon_i=1,...,W-1,\quad \psi_{\boldsymbol{\epsilon}\epsilon_i}=1$ when $\epsilon_i=W$ 

#### The lower bound of the difference between TSSB and TSSB-DW





### Equivalent form