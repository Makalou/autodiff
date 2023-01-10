# Automatic Differentiation

## Why AD programming matters?

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

---
## Todo List
- [ ] variadic arguments for ```gradient_at``` function
- [ ] optimization for reverse mode compute graph storage and traverse
- [ ] parallel computing

