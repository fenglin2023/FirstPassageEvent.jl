# FirstPassage.jl
This package provides a fast simulation algorithm for the first passage event of a subordinator, consisiting of a truncated tempered stable subordinator and a compound Poisson process, across a target function. 


This package includes many auxiliary functions, including the following functions. (The notation used below follows that of the original article published by the same authors of this package.)

* The following function produces a sample of $S_t$ under $\mathbb{P}_0$ (i.e. a stable increment or marginal):
```julia
rand_stable(α::Real, θ::Real, t::Real)
```

* The following function produces a sample of $S_t$ under $\mathbb{P}_q$ (i.e. a tempered stable increment or marginal):
```julia
rand_tempered_stable(α::Real, θ::Real, q::Real, t::Real)
```

* The following function produces a sample of $S_t|\\{S_t\le s\\}$ under $\mathbb{P}_0$:
```julia
rand_small_stable(α::Real, θ::Real, t::Real, s::Real)
```

* The following function produces a sample of $S_{t}|\\{S_{t}\le s\\}$ under $\mathbb{P}_{q}$:
```julia
rand_small_tempered_stable(α::Real, θ::Real, q::Real, t::Real, s::Real)
```

* The following function produces a sample of the undershoot $S_{\tau-}|\\{S_{\tau-} < S_{\tau}\\}$ under $\mathbb{P}_{0}$ where $\tau$ is the crossing time and $s = b(\tau)$ is the crossing level 
```julia
rand_undershoot_stable(α::Real, θ::Real, t::Real, s::Real)
```

* The following function produces a sample of the vector $(\tau, S_{\tau-}, S_{\tau}-S_{\tau-})$ under $\mathbb{P}_0$ where $\tau$ is the crossing time across the boundary $b$, `b` is the boundary function $b$, `Db` is the derivative of `b` and `B` is the inverse function of $t \mapsto t^{-1/\alpha}b(t)$: 
```julia 
rand_crossing_stable(α::Real, θ::Real, b::Function, Db::Function, B::Function)
```

* The following function produces a sample of the vector $(\tau, S_{\tau-}, S_{\tau}-S_{\tau-})|\\{t\le T\\}$ under $\mathbb{P}_0$ where $\tau$ is the crossing time across the boundary $b$, `T` is some positive number, `b` is the boundary function $b$, `Db` is the derivative of `b` and `B` is the inverse function of $t \mapsto t^{-1/\alpha}b(t)$:
```julia
rand_crossing_small_stable(α::Real, θ::Real, b::Function, Db::Function, B::Function, T::Real)
```

* The following function produces a sample of the vector $(\tau, S_{\tau-}, S_{\tau}-S_{\tau-})$ under $\mathbb{P}_q$ where $\tau$ is the crossing time across the boundary $b$ and $b:t\mapsto \min\\{a_0-a_1 t,r\\}$ is the boundary function:
```julia
rand_crossing_tempered_stable(α::Real, θ::Real, q::Real, a0::Real, a1::Real, r::Real)
```

* The following function produces the total mass $\Lambda$ of the compound Poisson component of the subordinator:
```julia
function Λ(α::Real,ϑ::Real,q::Real,r::Real,r0::Real) 
```

* The following function produces a sample of the jump $J$ of the compound Poisson component $Q$, $J$ has density proportional to $1{r0>x>r}e^(-qx)x^(-α-1)$ :
```julia
function randΛ(α::Real,q::Real,r0::Real,r::Real)
```

* The following function produces a sample of the vector $(\tau, Z_{\tau-}, Z_{\tau}-Z_{\tau-})$ where $\tau$ is the crossing time of subordinator $Z$ across the boundary $c$, $\rho$ is some number between $0$ and $1$:
```julia
function rand_crossing_subordinator(α::Real,ϑ::Real,q::Real,a0::Real,a1::Real,r::Real,ρ::Real,r0::Real,Λ0::Real,randΛ0::Function)
```
