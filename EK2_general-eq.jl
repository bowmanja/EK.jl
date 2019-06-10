###################################################################
#Description
###################################################################


###################################################################
##Packages
###################################################################
using LinearAlgebra, SpecialFunctions, NLsolve


###################################################################
##Parameters
###################################################################
σ = 3 #Elasticity of substitution, Source : Footnote 45, EK2002
θ = 8.28 #Technology dispersion parameter(Comparative Advantage), Source : Table VI, EK2002
β = 0.21 #Labor share in costs, Source : Table VIII, EK2002


###################################################################
#Generate data
#Wages, Trade costs, Technology shift parameter
###################################################################
N=43 #Number of countries

#Generate wages costs
ϵ_temp = randn(N)
w = 100*fill(1.0,N,1) + 100*ϵ_temp.^2

#Generate trade costs
ϵ_temp = randn(N,N)
D = fill(1.0,N,N) + ϵ_temp.^2
D=D-Diagonal(D)+I # Make diagonal entries equal 1 (No internal resistance)
#D = fill(1.0,N,N)

#Generate productivity shift parameter
ϵ_temp = randn(N)
T = 10*fill(1.0,N,1) + ϵ_temp.^2

###################################################################
#Model
###################################################################

############### Step 1 : Fixed point for prices ###################
function P_fixedpoint(P_telde,D,T,w,θ,β,σ,N)
    γ = (gamma((θ+1-σ)/θ))^(1/(1-σ)) #Gamma function
    Part1 = γ^(-θ)*T.*(w.^(-θ*β)) #create vector with product of gamma, wages^(θβ) and shift parameter
    Part2 = repeat(transpose(Part1),N) #create matrix which contain copies of above vector in each row
    G = (D.^(-θ)).*Part2 # point wise multiplication of distance^(-θ) and above matrix
    P_tilde = G*P_telde.^(1-β)
    return P_tilde
end
f(P_tilde) =  P_fixedpoint(P_tilde,D,T,w,θ,β,σ,N) #redefine function with single argument

init=ones(N,1) #initial value
sol = fixedpoint(f, init)
P = sol.zero.^(-1/θ)

############### Step 2 : Trade Shares ###################

γ = (gamma((θ+1-σ)/θ))^(1/(1-σ)) #Gamma function
Part1 = γ^(-θ)*T.*(w.^(-θ*β)).*(P.^(-θ*(1-β))) #create vector with product of gamma, wages^(θβ), P^(-θ(1-β)) and shift parameter
Part2 = repeat(transpose(Part1),N) #create matrix which contain copies of above vector in each row
Part3=transpose(repeat(transpose(P.^(-θ)),N)) #create matrix which contain copies of price vector in each row
Π = (D.^(-θ)).*Part2./Part3 # point wise multiplication of distance^(-θ) and above matrixs


#Correlation check
corr = cor(T,diag(Π))
println(corr)

############### Step 3 : Counterfactual ###################
# RTA between two countries : reduce distance between single pair and look at change in trade share for all countries


############### Step 4 : Graphics ###################
#Trade flows in circle
