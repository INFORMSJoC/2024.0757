[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Construction of Value Functions of Integer Programs with Finite Domain

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Construction of Value Functions of Integer Programs with Finite Domain](https://doi.org/10.1287/ijoc.2024.0757) by Y. Huang and J. Zhang. 

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

## Data
[Trapp et al. (2013)](https://doi.org/10.1287/opre.1120.1156) proposed 26 test instance classes for two-stage stochastic integer programs with random right-hand sides. 

For Algorithm 1, we evaluate its performance on the original instances constructed by Trapp et al. These instances are available under the folder `data/Nonnegative instances`.

For Algorithm 2, we generate modified instances derived from Trapp’s benchmarks by perturbing the constraint matrix. Specifically, we randomly select a subset of columns and, within those columns, reverse the signs of certain positive elements. 

More precisely:

- Each column in the constraint matrix is independently selected with probability $p$.  
- For each selected column, every positive element is then independently flipped in sign (i.e., made negative) with the same probability $p$.

The resulting instances are organized as follows:

- When $p = 10$%, the modified instances are stored in the folder `data/instances with 10% negative elements`.
- When $p = 5$%, the modified instances are located in the folder `data/instances with 5% negative elements`.

## Replicating
The algorithms are implemented in C++, compiled with Microsoft Visual Studio 2022, and executed on a desktop machine equipped with a 3.2 GHz single-core CPU, 64 GB of RAM, and CPLEX 12.9.0.

