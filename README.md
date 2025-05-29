<h1 align="center">Supplemental Material</h1>
<h2 aligne="center">Simple and accurate complete elliptic integrals for the full range of modulus</h2>
Teepais Chachiyo <teepanisc@nu.ac.th>, Department of Physics, Faculty of Science, Naresuan University, Phitsanulok 65000, Thailand.


## File List
1. **supplemental_material.ipynb**: a python notebook with additional proofs and figures as shown in the paper
2. **elliptic.py**: a python implementation of the formulae in this work
3. **elliptic.c**: a c language implementaion of the formulae in this work

<br>
The research article preprint >> <a href="https://arxiv.org/abs/2505.17159">https://arxiv.org/abs/2505.17159</a>

## Abstract
The complete elliptic integral of the first and second kind, K(k) and E(k), appear in a multitude of physics and engineering applications. Because there is no known closed-form, the exact values have to be computed numerically. Here, approximations for the integrals are proposed based on their asymptotic behaviors.  An inverse of K is also presented. As a result, the proposed K(k) and E(k) reproduce the exact analytical forms both in the zero and asymptotic limits, while in the mid-range of modulus maintain average error of 0.06\% and 0.01\% respectively.  The key finding  is the ability to compute the integrals with exceptional accuracy on both limits of elliptical conditions. An  accuracy of 1 in 1,000  should be sufficient for practical or prototyping engineering and architecture designs. The simplicity should facilitate discussions of advanced physics topics  in introductory physics classes, and enable broader collaborations among researchers from other fields of expertise. For example, the phase space of energy-conserving nonlinear pendulum using only elementary functions is discussed. The proposed inverse of K is shown to be Never Failing Newton Initialization and is an important step for the computation of the exact inverse. An algorithm based on Arithmetic-Geometric Mean for computing exact integrals and their derivatives are also presented. It should be useful in a platform that special functions are not accessible such as web-based and firmware developments.

## Citation

In the meantime, if you use any part of this repository please cite the following preprint:

```
@article{Chachiyo:2025saa,
    author = "Teepanis Chachiyo",
    title = "{Simple and accurate complete elliptic integrals for the full range of modulus}",
    eprint = "2505.17159",
    archivePrefix = "arXiv",
    primaryClass = "physics.class-ph",
    month = "4",
    year = "2025",
    note = {arXiv:2505.17159},
}
