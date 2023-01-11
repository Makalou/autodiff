# Automatic Differentiation

## Why AD programming matters?
[Automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)(AD), or autodiff for short, can evaluate the 
derivative of a function specified by a computer program automatically. Which means, once you write down your function in C++ code, 
you don't have to suffer from the intricate and error-prone derivative deduction on the paper, then code the result manually. 
All you need to do is passing the function to autodiff library, and get its derivative or gradient whenever you want :) !
### Optimization & Gradient Descent Method
Once you write down a function, looking for the parameters that make the maximum or minimum value, you are involved in an 
[optimization problem](https://en.wikipedia.org/wiki/Optimization_problem). Generally speaking, an optimization problem is the problem of finding the 
best solution from all feasible solutions. For example, blabla...
There is one thing worth noting: depending on the whether the variables are discrete or continuous, optimization can be divided into two categories:
[discrete optimization](https://en.wikipedia.org/wiki/Discrete_optimization) and [continuous optimization](https://en.wikipedia.org/wiki/Continuous_optimization). 
In discrete optimization, both the inputs and outputs are discrete. To solve such problem, combinatorial mathematics and integer programming can come to the stage.
The other form of the optimization, the continuous one, expects the inputs and outputs changing continuously and smoothly. One of the most famous algorithms to solve
such problems is [gradient descent](https://en.wikipedia.org/wiki/Gradient_descent). 
### Other Applications

## Automatic Differentiation Mode

### Forward mode & Dual number

### Reverse mode & Compute graph
### Forward or Reverse?
- Do forward mode and reverse mode have any numeric accuracy difference?
  - No.
- Which mode should I choose when the target function has both inputs and outputs in really high dimension?
  - Under this case unfortunately neither **forward** nor **reverse** mode could be expected to execute efficiently. 
    The computing graph is complex so neither forward nor backward traverse can help.

## High Performance Computing
### Multi-threading
### SIMD
### GPU Acceleration

---
## Todo List
- [ ] variadic arguments for ```gradient_at``` function
- [ ] optimization for reverse mode compute graph storage and traverse
- [ ] parallel computing

