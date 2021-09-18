Hey,

A Complementary software for the paper: "The generalized method of moments for multi-reference alignment" by Asaf Abas, Tamir Bendory,	and~Nir Sharon.

This code is a part of on going to improve the reconstruction results in cryo-EM.
I have implemented the generlized method of moment (GMM) for a quick intial reconstruciton of the volume.

However, this code is not stand alone. It is an expansion of the work in: https://web.math.princeton.edu/~amits/publications/1907.05377.pdf    (Method of moments for 3-D single particle ab initio modeling with non-uniform distribution of viewing angles)

In order to use the code: 
1. Download the full framework from: https://github.com/nirsharon/nonuniformMoM
2. Install/download all the third party software as writen in https://github.com/nirsharon/nonuniformMoM
3. Add the GMM folder the main directory of full frameworks
4. Add the GMM folder to the new path in matlab
5. Run the GMM code via one of the file in the new GMM folder, sub-folder: TestCase

Best,
Asaf


The reconstruciton results in the paper are for:
- seed: 10
-grid size: 23
- number of observaitons: 20,000
- sigma = 0.05
- Exp numeber :48

