# General

This repo provides a Matlab interface to the [HiGHS](https://highs.dev/) solver. Two drivers are generated:
- `highs_lp` a miLP interface
- `highs_qp` a QP interface

## The linear programming driver

The LP interface solves the [problem formulation](https://ergo-code.github.io/HiGHS/dev/#Specification). We pass

```
[x,fVal,info] = highs_lp(c,A,[L,U])
[x,fVal,info] = highs_lp(c,A,[L,U],[l,u])
[x,fVal,info] = highs_lp(c,A,[L,U],[l,u],isIntegral)
```

where `c`, `A`, `L`, `U`, `l`, `u` have the meaning specified in [the problem formulation](https://ergo-code.github.io/HiGHS/dev/#Specification) and `isIntegral` is a `logical` array where `true` implies that an integer solution should be obtained.
Missing bounds are passed by specifying them as `inf` or `-inf` accordingly.

## The quadratic programming driver

The QP driver does not admit solving for mixed integer problems see the [HiGHS documetation](https://ergo-code.github.io/HiGHS/dev/#Specification). The remainder of the interface is analogous to the LP interface

```
[x,fVal,info] = highs_qp(Q,c,A,[L,U])
[x,fVal,info] = highs_qp(Q,c,A,[L,U],[l,u])
```

where the only novelty is that the Hessian is passed as a triangular matrix.

## The return structure

The `info` structure is contains the dual variables for both the constraint rows `A` as well as the bounds on the decision variable (denoted column within HiGHS). Furthermore, additional information on the status of the respective constraint is returned.


## Building

The project is entirely based on `cmake` and should build straight out-of-the-box, given that a current version of cmake is installed. 
Unit tests, as well as presets for the main platforms are available.

