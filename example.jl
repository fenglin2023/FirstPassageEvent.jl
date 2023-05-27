#OU process
#this is the code of the example in section 4.2
using Plots, Random, SpecialFunctions, Distributions, LaTeXStrings

Random.seed!(2023)

rand_crossing_tempered_stable(.55, 1., 3, 4, 0, 2)

α = 0.65
ϑ = 1/(-gamma(-α))
q = 1.
r = 1

n = Int(1e3)
xarr = convert(Array{Float64},0:.05:4)
yarr = convert(Array{Float64},0:.05:4)
T, U = zeros(n), zeros(n)

u = zeros(length(xarr),length(yarr))
gmma = [1 1; 0 1]#gmma is the volatility of the SDE
function myrandΛ()# Q has levy density 1{s>1}s^(-5)
    return rand()^(-1/4)
end
for k in (1:1:n)
    if k%10 == 0
        GC.gc()
        println("Iteration for $k")
    end
    t = 5# t is the  boundary 
    if U[k] < t
        (auxT,auxU,auxV) = rand_crossing_subordinator(α, ϑ, q, t-U[k], 0, 1,1/2,1/4,myrandΛ) 
        T[k] += auxT
        U[k] += auxU + auxV
    end
    for i in eachindex(xarr)
        for j in eachindex(yarr)
            dd = MvNormal([0;0],gmma*transpose(gmma).*((1-exp(-2*T[k])/2)))
            x = exp(-[1 0; 0 1].*T[k])*[i;j]+rand(dd,1)
            u[i,j] += x[1,1]+x[2,1]^2
        end
    end
end
u ./= n

surface(xarr, yarr, transpose(u), xlabel = L"$x_1$", ylabel = L"$x_2$", zlabel = L"$u$", title = L"Solution $u(t,x)$ of the FPDE")