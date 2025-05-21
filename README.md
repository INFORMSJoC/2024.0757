[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0757

https://doi.org/10.1287/ijoc.2024.0757.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Huang2025,
  author =        {Yifei Huang, Junlong Zhang},
  publisher =     {INFORMS Journal on Computing},
  title =         {Construction of Value Functions of Integer Programs with Finite Domain},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0757.cd},
  url =           {https://github.com/INFORMSJoC/2024.0757},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0757},
}  
```

## Description

Value functions play a central role in integer programming duality, and they are also used to develop solution methods for stochastic integer programs, bilevel integer programs and robust optimization problems. In this paper, we propose a column-by-column algorithm for constructing the value functions of integer programs with finite domains over the set of level-set minimal vectors. The proposed algorithm starts with the first column and sequentially adds the rest of the columns one by one. Each time a column is added, a new set of level-set minimal vectors is generated based on the previous set, and the optimal objective values over the level-set minimal vectors are also computed. The advantage of the proposed algorithm is that no integer program needs to be solved in the algorithm for instances with nonnegative constraint matrices. Computational results on benchmark instances show that the proposed algorithm can achieve a speedup of up to three orders of magnitude compared with a state-of-the-art algorithm. We also extend the proposed algorithm to build value functions of integer programs with negative elements in the constraint matrix.
