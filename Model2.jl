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
θ = 8.28 #Technology dispersion parameter(Comparative Advantage), Source : Table VI, EK2002
β = 0.21 #Labor share in costs, Source : Table VIII, EK2002


###################################################################
#Generate data
#Wages, Trade costs, Technology shift parameter
###################################################################
N=10 #Number of countries

#Generate wages costs
ϵ_temp = randn(N)
w = 1000*fill(1.0,N,1) + ϵ_temp.^2

#Generate trade costs
ϵ_temp = randn(N,N)
D = fill(1.0,N,N) + ϵ_temp.^2
D=D-Diagonal(D)+I # Make diagonal entries equal 1 (No internal resistance)
#D = fill(1.0,N,N)

#Generate productivity shift parameter
ϵ_temp = randn(N)
T = 1000*fill(1.0,N,1) + ϵ_temp.^2

###################################################################
#Model
###################################################################

#price function
function prices(P_telda,D,T,w,θ,β,σ,N)
    #Gamma function
    γ = (gamma((θ+1-σ)/θ))^(1/(1-σ))

    Part1 = γ^(-θ)*T.*(w.^(-θ*β)) #create vector with product of gamma, wages^(θβ) and shift parameter
    Part2 = repeat(transpose(T),N) #create matrix which contain copies of above vector in each row
    G = (D.^(-θ)).*Part2 # point wise multiplication of distance^(-θ) and above matrix
    error_vec = P_telda - G*P_telda # calculate error in current equation
    error = transpose(error_vec)*error_vec
    return error
end


f = P_tilde->prices(P_tilde,D,T,w,θ,β,σ,N) #redefine function with single argument


#=
find how to minimize f function
=#

b=f(fill(1.0,N,1))


using NLsolve

function g!(a,P_tilde)
    a = f(P_tilde)
    println(a[1])
end
x_0 = fill(4.0,N,1)
a=nlsolve(g!,x_0,autodiff = :forward)




using Optim
a = optimize(f,fill(1.0,N,1))
