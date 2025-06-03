> [!NOTE] 
> For proper equation rendering, please view this documentation in day mode instead of night mode. 

This repository provides the MATLAB implementation of the MCMC algorithm presented in the following paper: 

[Takashi Furuya](https://kendb.doshisha.ac.jp/profile/en.77bfc1f47b9eacdc.html), [Pu-Zhao Kow](https://puzhaokow1993.github.io/homepage/) and [Jenn-Nan Wang](http://www.math.ntu.edu.tw/~jnwang/), *Consistency of the Bayes method for the inverse scattering problem with randomly truncated Gaussian priors*, May 14, 2025, [preprint](https://www.math.ntu.edu.tw/~jnwang/pub/resources/random_truncation(0515).pdf)

In this work, we apply the Bayes approach to study the inverse medium scattering problem. Let ![n\ge0](https://latex.codecogs.com/png.image?\dpi{110}n\ge0) and ![1-n](https://latex.codecogs.com/png.image?\dpi{110}1-n) be a compactly supported function in ![\mathbb{R}^{3}](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^{3}) with ![{\rm%20supp}%20(1-n)\subset%20D](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20supp}%20(1-n)\subset%20D), where ![D](https://latex.codecogs.com/png.image?\dpi{110}D) is an open bounded smooth domain, and having suitable regularity. Let ![u_n%20=u^{\rm%20inc}+u_n^{\rm%20sca}](https://latex.codecogs.com/png.image?\dpi{110}u_n%20=u^{\rm%20inc}+u_n^{\rm%20sca}) satisfy 
<div align="center">
  
![\Delta%20u_n+k^2%20nu_n=0](https://latex.codecogs.com/png.image?\dpi{110}\Delta%20u_n+k^2%20nu_n=0) in ![\mathbb{R}^3](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^3)
</div>

and the scattered field ![u_n^{\rm%20sca}](https://latex.codecogs.com/png.image?\dpi{110}u_n^{\rm%20sca}) satisfies the Sommerfeld radiation condition 
<div align="center">
  
![\lim_{|x|\rightarrow\infty}|x|\left(\frac{\partial%20u_n^{\rm%20sca}}{\partial|x|}-\mathbf{i}ku_n^{\rm%20sca}\right)=0](https://latex.codecogs.com/png.image?\dpi{110}\lim_{|x|\rightarrow\infty}|x|\left(\frac{\partial%20u_n^{\rm%20sca}}{\partial|x|}-\mathbf{i}ku_n^{\rm%20sca}\right)=0).  
</div>

Assume that ![u^{\rm%20inc}](https://latex.codecogs.com/png.image?\dpi{110}u^{\rm%20inc}) is the plane incident field, i.e. ![u^{\rm%20inc}(\theta',\theta)=e^{\mathbf{i}k\theta'\cdot\theta}](https://latex.codecogs.com/png.image?\dpi{110}u^{\rm%20inc}(\theta',\theta)=e^{\mathbf{i}k\theta'\cdot\theta}) for all ![\theta',\theta\in\mathbb{S}^2](https://latex.codecogs.com/png.image?\dpi{110}\theta',\theta\in\mathbb{S}^2). Then the scattered field ![u_n^{\rm%20sca}](https://latex.codecogs.com/png.image?\dpi{110}u_n^{\rm%20sca}) posses the asymptotic behavior 
<div align="center">
  
![u_n^{\rm%20sca}(x,\theta)=\frac{e^{\mathbf{i}k|x|}}{|x|}u_n^\infty(\theta',\theta)+o(r^{-1})](https://latex.codecogs.com/png.image?\dpi{110}u_n^{\rm%20sca}(x,\theta)=\frac{e^{\mathbf{i}k|x|}}{|x|}u_n^\infty(\theta',\theta)+o(r^{-1}))   as   ![|x|\rightarrow\infty](https://latex.codecogs.com/png.image?\dpi{110}|x|\rightarrow\infty),  
</div>

where ![\theta'=x/|x|](https://latex.codecogs.com/png.image?\dpi{110}\theta'=x/|x|). The inverse scattering problem is to determine the medium perturbation ![1-n](https://latex.codecogs.com/png.image?\dpi{110}1-n) from the knowledge of the scattering amplitude ![u_n^{\infty}(\theta',\theta)=e^{\mathbf{i}k\theta'\cdot\theta}](https://latex.codecogs.com/png.image?\dpi{110}u_n^{\infty}(\theta',\theta)) for all ![\theta',\theta\in\mathbb{S}^2](https://latex.codecogs.com/png.image?\dpi{110}\theta',\theta\in\mathbb{S}^2). We adopt a Bayesian approach to this problem, providing not only rigorous theoretical justifications but also supporting numerical simulations.

# Algorithm # 

**Require:** an initial guess ![F^{(0)}](https://latex.codecogs.com/png.image?\dpi{110}F^{(0)}) and a resolution parameter ![J_{(0)}](https://latex.codecogs.com/png.image?\dpi{110}J_{(0)})

1. Set ![{\rm%20status}(0)={\rm%20true}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(0)={\rm%20true})
2. Creating an empty sequence ![\{{\rm%20status}(\tau)\}_{\tau=1}^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}\{{\rm%20status}(\tau)\}_{\tau=1}^{\infty})

# Results # 

TBA

[comment]: <> (https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)
