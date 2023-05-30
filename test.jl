#this test compare the empirical cdfs of the first passage time of subordiantor
#with Levy measure 1{0<x<1}e^(-x)x^(-α-1)dx+[1{0<x<1}(1-e^(-x))x^(-α-1)dx+1{x>1}x^{-α-1}dx],
#using our algorithm, and using the law of (c/S(1))^α.
using Plots, SpecialFunctions, StatsBase,  Random
Random.seed!(2022)

α = .6
ϑ = 1.
b = 1.
q = 1.
a0 = 3.
r = 1
Λ1 = 2.20171#\int_0^1 (1-e^(-qx))x^(-α-1)dx
Λ2 = 1.66667#\int_0^1 x^(-α-1)dx

function my_rand()#sample from the jump law vartheta (1-e^{-qx})x^{-alpha-1}dx
    if rand()<Λ1/(Λ1+Λ2)
        X = rand()^(1/(1-α))
        while rand()>(1-exp(-q*X))/(q*X)
            X = rand()^(1/(1-α))
        end
    else
        X = rand()^(-1/α)
    end
    return X
end

n = 10000

FPE =[rand_crossing_subordinator(α,ϑ,q,a0,0,r,1/2,Λ1+Λ2,my_rand) for i = 1:n]
passage_time = map(x -> x[1], FPE)
XT = collect(passage_time)
FT = [ (a0/rand_stable(α,ϑ*gamma(1-α)/α,1))^α for i = 1:n] #the time has the same law of (a0/S(1))^α
# Empirical CDFs
ECDF_T = ecdf(XT)
ECDF_t = ecdf(FT)
EFT(x) = ECDF_T(x)
EFt(x) = ECDF_t(x)
plot([EFT, EFt], 0.01, 5)

# scaled difference between ECDFs
BB(x) = sqrt(n/2)*(ECDF_T(x)-ECDF_t(x))
xarr = 0:0.01:3
plot(ECDF_t(xarr),BB(xarr),label= "scaled difference between ECDFs", ylim=(-1,1))
