# HPhiJulia
Julia wrapper for the HPhi.
This just calls HPhi from Julia.
Now this is very limited. 
I confirmed that this can work on Ubuntsu and Mac OS 10.13. 

# HPhi 
Quantum Lattice Model Simulator Package
https://github.com/issp-center-dev/HPhi

HPhi is installed in ".julia/packages/HPhiJulia".


# Install 


```julia
add https://github.com/cometscome/HPhiJulia
```

#Examples

## Heisenberg on a 2D square lattice

```julia
using HPhiJulia
L = 4
W = 4
J = 1
results = HPhiJulia.HPhi("Spin","square",W,L,J=J,mpinum=1)
println(results["Energy"])
```

## Hubbard model on a 2D square lattice

```julia
using HPhiJulia
L = 2
W = 2
U = 8
t = 1
nelec = 4
results =  HPhiJulia.HPhi("Hubbard","square",W,L,U=U,t=t,nelec=nelec)
println(results["Energy"])
```
