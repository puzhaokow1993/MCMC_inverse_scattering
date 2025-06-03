> [!NOTE] 
> For proper equation rendering, please view this documentation in day mode instead of night mode. 

This repository provides the MATLAB implementation of the MCMC algorithm presented in the following paper: 

[Takashi Furuya](https://kendb.doshisha.ac.jp/profile/en.77bfc1f47b9eacdc.html), [Pu-Zhao Kow](https://puzhaokow1993.github.io/homepage/) and [Jenn-Nan Wang](http://www.math.ntu.edu.tw/~jnwang/), *Consistency of the Bayes method for the inverse scattering problem with randomly truncated Gaussian priors*, May 14, 2025, [preprint](https://www.math.ntu.edu.tw/~jnwang/pub/resources/random_truncation(0515).pdf)

We will not going to explain all notations here, please refer to our paper for more details. 

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

1. Set ![{\rm%20status}(0)={\rm%20accept}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(0)={\rm%20accept})
2. Creating an empty sequence ![({\rm%20status}(\tau))_{\tau=1}^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}({\rm%20status}(\tau))_{\tau=1}^{\infty})
3. **for ![\tau=0,1,2,\cdots](https://latex.codecogs.com/png.image?\dpi{110}\tau=0,1,2,\cdots) do**
4. $~~~~$ Generate ![F(x,y)=\sum_{r,s\in\mathbb{Z}}F_{rs}\chi_{(0,1)\times(0,1)}(2^{J_0}x-r,2^{J_0}y-s)](https://latex.codecogs.com/png.image?\dpi{110}F(x,y)=\sum_{r,s\in\mathbb{Z}}F_{rs}\chi_{(0,1)\times(0,1)}(2^{J_0}x-r,2^{J_0}y-s)) with randomly chosen ![F_{rs}\sim\mathcal{N}(0,1)](https://latex.codecogs.com/png.image?\dpi{110}F_{rs}\sim\mathcal{N}(0,1))
5. $~~~~$ Propose ![F^{(\tau+1)}=\sqrt{(1+\beta^2)}F^{(\tau)}+\beta%20F](https://latex.codecogs.com/png.image?\dpi{110}F^{(\tau+1)}=\sqrt{(1+\beta^2)}F^{(\tau)}+\beta%20F)
6. $~~~~$ **if ![{\rm%20status}(\tau)={\rm%20accept}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau)={\rm%20accept}) then**
7. $~~~~~~~~$ ![\ell_{\rm%20current}=\ell^{(N)}(F^{(\tau)})](https://latex.codecogs.com/png.image?\dpi{110}\ell_{\rm%20current}=\ell^{(N)}(F^{(\tau)})), where ![\ell^{(N)}](https://latex.codecogs.com/png.image?\dpi{110}\ell^{(N)}) is the log-likelihood function given in equation (2.8) of our paper.
8. $~~~~$ **end if**
9. $~~~~$ Compute ![u_*=\ell^{(N)}(F^{(\tau+1)})-\ell_{\rm%20current}](https://latex.codecogs.com/png.image?\dpi{110}u_*=\ell^{(N)}(F^{(\tau+1)})-\ell_{\rm%20current})
10. $~~~~$ Generate a uniform random number ![u\in[0,1]](https://latex.codecogs.com/png.image?\dpi{110}u\in[0,1])
11. $~~~~$ **if ![log(u)<u_*](https://latex.codecogs.com/png.image?\dpi{110}log(u)<u_*) then**
12. $~~~~~~~~$ Set ![{\rm%20status}(\tau+1)={\rm%20accept}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau+1)={\rm%20accept})
13. $~~~~$ **else**
14. $~~~~~~~~$ Set ![F^{(\tau+1)}=F^{(\tau)}](https://latex.codecogs.com/png.image?\dpi{110}F^{(\tau+1)}=F^{(\tau)}) and ![{\rm%20status}(\tau+1)={\rm%20reject}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau+1)={\rm%20reject})
15. $~~~~$ **end if**
16. **end for**

**Return:** ![(F^{(\tau)})_{\tau=1}^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}(F^{(\tau)})_{\tau=1}^{\infty}) by removing all entries ![F^{(\tau)}](https://latex.codecogs.com/png.image?\dpi{110}F^{(\tau)}) each corresponds to ![{\rm%20status}(\tau)={\rm%20reject}](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20status}(\tau)={\rm%20reject})

> [!IMPORTANT] 
> The if condition ![log(u)<u_*](https://latex.codecogs.com/png.image?\dpi{110}log(u)<u_*) and the log-likelihood function ![\ell^{(N)}](https://latex.codecogs.com/png.image?\dpi{110}\ell^{(N)}) cannot be replaced by ![u<\exp(u_*)](https://latex.codecogs.com/png.image?\dpi{110}u<\exp(u_*)) and the likelihood function, respectively, due to limited floating point, despite they are mathematical equivalent. 

> [!TIP] 
> The numerical computation can be significantly speeded up by employing the parallel computing, using [MATLAB Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html). Using `parfor` rather than `for` also speed up the computations. 

# Results # 

The perturbation ![1-n^{(\rm%20true)}](https://latex.codecogs.com/png.image?\dpi{110}1-n^{(\rm%20true)}) is approximated by ![1-\Phi(\overline{F}_{\rm%20burn}^{(T)})](https://latex.codecogs.com/png.image?\dpi{110}1-\Phi(\overline{F}_{\rm%20burn}^{(T)})) for a sufficiently large integer ![T](https://latex.codecogs.com/png.image?\dpi{110}T), where ![\overline{F}_{\rm%20burn}^{(T)}](https://latex.codecogs.com/png.image?\dpi{110}\overline{F}_{\rm%20burn}^{(T)}) denotes the 'burn-in' sample mean, defined as 
<div align="center">
  
![\overline{F}_{\rm%20burn}^{(T)}=\frac{1}{\lfloor%20T/2\rfloor}\sum_{\tau=\lfloor%20T/2\rfloor+1}^{T}F^{(\tau)}](https://latex.codecogs.com/png.image?\dpi{110}\overline{F}_{\rm%20burn}^{(T)}=\frac{1}{\lfloor%20T/2\rfloor}\sum_{\tau=\lfloor%20T/2\rfloor+1}^{T}F^{(\tau)})
</div>

The calculations were performed on an 8-core workstation, allowing a parallel solution of 8 equations simultaneously. The PDE (used to compute the log-likelihood function) was discretized using the finite element method with a mesh size of at most 0.3. 

## Result 1 ##

In order to reconstruct 
<div align="center">
  
![n^{(\rm%20true)}=1+\frac{1}{10}\chi_{[-1/2,1/2]\times[-1/2,1/2]}](https://latex.codecogs.com/png.image?\dpi{110}n^{(\rm%20true)}=1+\frac{1}{10}\chi_{[-1/2,1/2]\times[-1/2,1/2]}), 
</div>

we choose the wave number ![k=5](https://latex.codecogs.com/png.image?\dpi{110}k=5), noise level ![\sigma=0.1](https://latex.codecogs.com/png.image?\dpi{110}\sigma=0.1), initial guess ![F^{(0)}\equiv0](https://latex.codecogs.com/png.image?\dpi{110}F^{(0)}\equiv0), sample size ![N=560](https://latex.codecogs.com/png.image?\dpi{110}N=560), learning rate ![\beta=0.05](https://latex.codecogs.com/png.image?\dpi{110}\beta=0.05), and resolution parameter ![J_0=1](https://latex.codecogs.com/png.image?\dpi{110}J_0=1). A total of 30000 iterations were carried out, requiring approximately a duration of 18 days and 3 hours (specifically, 1568082 seconds) to compute. 

The progression of the log-likelihood ![\ell_{\rm%20current}](https://latex.codecogs.com/png.image?\dpi{110}\ell_{\rm%20current}) over iterations is shown in figure below. 

<div align="center">
<img src="square_log_likelihood.png" alt="Our results" width="45%" style="margin-right: 10px;" /> 
</div>

The plot of 'burn-in' sample mean ![\overline{F}_{\rm%20burn}^{(T)}](https://latex.codecogs.com/png.image?\dpi{110}\overline{F}_{\rm%20burn}^{(T)}) and the ground truth ![n^{(\rm%20true)}](https://latex.codecogs.com/png.image?\dpi{110}n^{(\rm%20true)}) on FEM mesh is shown in figure below. 

<div align="center">
<img src="square_plot1_from_top.png" alt="Our results" width="45%" style="margin-right: 10px;" /> <img src="square_true_n_from_top.png" alt="Our results" width="45%" style="margin-right: 10px;" /> 
</div>

## Result 2 ##

In order to reconstruct 
<div align="center">
  
![n^{(\rm%20true)}=1+\frac{1}{10}\chi_{x<1/2,y<1/2,x+y>-1/4}](https://latex.codecogs.com/png.image?\dpi{110}n^{(\rm%20true)}=1+\frac{1}{10}\chi_{x<1/2,y<1/2,x+y>-1/4}), 
</div>

we choose the wave number ![k=5](https://latex.codecogs.com/png.image?\dpi{110}k=5), noise level ![\sigma=0.1](https://latex.codecogs.com/png.image?\dpi{110}\sigma=0.1), initial guess ![F^{(0)}\equiv0](https://latex.codecogs.com/png.image?\dpi{110}F^{(0)}\equiv0), sample size ![N=560](https://latex.codecogs.com/png.image?\dpi{110}N=560), learning rate ![\beta=0.05](https://latex.codecogs.com/png.image?\dpi{110}\beta=0.05), and resolution parameter ![J_0=1](https://latex.codecogs.com/png.image?\dpi{110}J_0=1). Based on some experiences in [Result 1](#result-1), this time only a total of 10000 iterations were carried out, requiring approximately a duration of 5 days and 8 hours (specifically, 459934 seconds) to compute. 

The progression of the log-likelihood ![\ell_{\rm%20current}](https://latex.codecogs.com/png.image?\dpi{110}\ell_{\rm%20current}) over iterations is shown in figure below. 

<div align="center">
<img src="triangle_log_likelihood.png" alt="Our results" width="45%" style="margin-right: 10px;" /> 
</div>

The plot of 'burn-in' sample mean ![\overline{F}_{\rm%20burn}^{(T)}](https://latex.codecogs.com/png.image?\dpi{110}\overline{F}_{\rm%20burn}^{(T)}) and the ground truth ![n^{(\rm%20true)}](https://latex.codecogs.com/png.image?\dpi{110}n^{(\rm%20true)}) on FEM mesh is shown in figure below. 

<div align="center">
<img src="triangle_plot1_from_top.png" alt="Our results" width="45%" style="margin-right: 10px;" /> <img src="triangle_true_n_from_top.png" alt="Our results" width="45%" style="margin-right: 10px;" /> 
</div>


[comment]: <> (https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)
