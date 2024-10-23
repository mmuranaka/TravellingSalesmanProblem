# Travelling Salesman Problem

## 1. Overview

This project involved finding optimized solutions to the Travelling Salesman Problem. The program receives a list of places with coordinates and attempts to solve a path to find the least expensive sollution. Two different algorithms, genetic and simulated annealing, were built.

## 2. Algorithms

**Genetic Algorithm**

The genetic algorithm begins with an initial population of random solutions and aims to evolve between generations into a more optimized solution through selective breeding of fit solutions. Of the initial population 20% of the population (those with that are most optimized) are elitist and does not have to go through the selection process. Of that 20%, half of the solutions will be untouched while the other half will undergo mutation. The population to undergo mutation will experience a switch between two random places along the salesman's path. The fitness function to determine parent likeliness is the total weight of the salesman's path while the selection mechanism involves randomly selecting a parent where more fit parents are most likely to be picked. The breeding between parents involved a two-point crossover, a child can consist of the beginning and end of parent 1 and the middle of parent 2.

Key Ideas:
* *pool* is the initial population size, I have set the defualt population size to 10
* Fitness Function the total weight of a solution
* Elitism is experienced by 20% of the population
  * half of elitists move on to the next generation untouched
  * the other half undergo mutation
* Mutation involves switching two random places in the path
  * Probability is given by the user
* Selection Mechanism is roulette wheel where most fit solutions are more likely to be randomly chosen for breeding
* Crossover Mechanism is two-point crossover where selected sections are switched between the parents
* Termination Condition is decided by the number of iterations specified by the user


**Simulated Annealing**

The simulated annealing algorithm begins with an intitial solution chosen at random and an initial temperature given to the algorithm by the user. To find an optimized solution to the Travelling Salesman problem, the algorithm will perform a two edge swap, similar to the mutation mechanism described for the genetic algorithm. If the fitness of this new solution is greater than that of the current solution, this new solution becomes our current solution. Otherwise, there is a probability that this new solution will replace the current solution where the probability is based on the objective function. Each iteration of this process cools the initial temperature and the algorithm will generate new solutions until the temperature drops below zero. The current solution following this termination condition is the algorithm's optimized solution.

Key Ideas:
* Initial Temperature is given by the user
* Fitness Function is the total cost of the path
* Objective Function is $P = e^{ \Delta Fitness / \ T_i}$
  * probability that the new, worse solution will replace the current solution
* Termination Condition is when the temperature is less than or equal to 0
* Cooling Schedule is exponential 
  * $T_i = T_{\text{i-1}} \times e^{-i \times \alpha}$
  * *i* is the iterations
  * *alpha* is the cooling parameter given by the user
  
## 3. Usage

The program can be executed with the line:  
`python tester.py DATA.CSV ALGO P1 P2`

- `tester.py` is the name of the python program
- `DATA.CSV` is the name of the csv file containing the vertices to be solved
- `ALGO` specifies the algorithm to be used
  - `1` specifies simulated annealing
  - `2` specifies genetic algorithm
- `P1` is a parameter specifying:
  - initial temperature for simulated annealing
  - number of iterations for genetic algorithm
- `P2` is a parameter specifying:
  - cooling parameter for simulated annealing
  - mutation probability for genetic algorithm
