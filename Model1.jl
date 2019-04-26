###################################################################
#=
Description :
=#
###################################################################


###################################################################
##Packages
###################################################################
using LinearAlgebra, SpecialFunctions


###################################################################
##Parameters
###################################################################
σ = 3 #Elasticity of substitution, Source : Footnote 45, EK2002
θ = 8.28 #Technology dispersion parameter, Source : Table VI, EK2002


###################################################################
#Generate data
#Unit costs, Trade costs, Technology shift parameter
###################################################################
N=10 #Number of countries

#Generate unit costs
ϵ_temp = randn(N)
c = 100*fill(1.0,N,1) + ϵ_temp.^2
#Generate trade costs
ϵ_temp = randn(N,N)
#D = fill(1.0,N,N) + ϵ_temp.^2
#D=D-Diagonal(D)+I # Make diagonal entries equal 1 (No internal resistance)
D = fill(1.0,N,N)

#Generate productivity shift parameter
ϵ_temp = randn(N)
T = 1000*fill(1.0,N,1) + ϵ_temp.^2

###################################################################
#Model
###################################################################

#Multilateral Resistance
Part1 = T.*c.^(-θ)
Φ = (D.^(-θ))*Part1

#Gamma function
γ = (gamma((θ+1-σ)/θ))^(1/(1-σ))

#Price Index
P = γ*Φ.^(-1/θ)

#Trade shares
Π = Diagonal((Φ.^(-1))*fill(1.0,1,N))*(D.^(-θ))*Diagonal(Part1*fill(1.0,1,N))




covT =(1/N)*diag(Π)'T - ((1/N)*(transpose(diag(Π))*fill(1.0,N,1)))*((1/N)*(transpose(T)*fill(1.0,N,1)))
