# Automatic Differentiation

## Why AD programming matters?
[Automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)(AD), or autodiff for short, can evaluate the derivative of a function specified by a computer program
automatically. Which means, once you write down your function in C++ code, you don't have to suffer from the intricate and error-prone
derivative deduction on the paper, then code the result manually. All you need to do is passing the function to autodiff library, and 
get its derivative or gradient whenever you want :) !
### Optimization & Gradient Descent Method

### Other Applications

## Automatic Differentiation Mode

### Forward mode & Dual number

### Reverse mode & Compute graph
### Forward or Reverse?
- Do foward mode and reverse mode have any numeric accuracy difference?
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

