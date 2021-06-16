################# Real Business Cycle - based on the model by Brock and Mirmam (1972)
################# Jo?o Costa-Filho
################# sites.google.com/site/joaoricardocostafilho
################# Twitter: @costafilhojoao


#### Parameters of the model ####

alpha = 0.4   # capital share in the production function
delta = 0.1   # capital depreciation rate
gamma = 2.0   # coefficient of relative risk aversion
beta  = 0.96  # time discount variable
nk    = 1000  # number of k grid points
nz    = 7     # number of possible states for the economy
rho   = 0.9   # persistency of the shock
std   = 0.05  # standard deviation of the shocks
m     = 1     # the number of standard deviations to approximate
simul = 30    # number of simulation periods


#### Discretizing the z grid #### 
# we can approximate the AR(1) process of z by a markov chain using the Tauchen (1986)

library(Rtauchen)

zprob <- Rtauchen( nz, std, rho, m )        # probability transition matrix
z     <- exp ( Tgrid( nz, std, rho, m ) )   

#### Steady State ####

ky_ss = ( beta * alpha ) / ( 1 - beta * ( 1 - delta ) )
y_ss  =  ky_ss^( alpha / ( 1 - alpha ) )  
k_ss  = ky_ss * y_ss       
c_ss  = y_ss - delta * k_ss
z_ss  = 1                                 # the steady state value for z is normalized to unity

#### Discretizing the k grid #### 

kmin = 0.5 * k_ss
kmax = 1.5 * k_ss

k = matrix( seq( from = kmin, to = kmax-(kmax-kmin)/nk,
                 by = (kmax-kmin)/nk),
            nrow = 1, ncol = nk )


#### Functional forms ####


utility <- function( c, gamma ){
  
    "
    Calculates the utility given the level of consumption
    and the relative risk aversion coeficient

    Parameters
    ----------
        c: the level of consumption
    gamma: coefficient associated with the relative risk aversion
    "

    c^( 1 - gamma ) / ( 1 - gamma )
      
}

#### Initial values for the value function ####

u = utility(c_ss, gamma)                               # steady state utility
V = matrix(1, nrow = nz, ncol = nk) * u / ( 1 - beta ) # infinite sum of discounted utilities at constant steady state consumption

v     = matrix(0, nrow = nz, ncol = nk)                # Initializes value function
g     = matrix(0, nrow = nz, ncol = nk)                # Initializes policy function
newV  = matrix(0, nrow = nz, ncol = nk)                # Initializes first step iteration


#### Value function iteration --- large nk method (grid search) ####

# Tolerance levels

tol    = 10^-4
maxiter = 300


iter = 0
d    = 1

while (d > tol && iter < maxiter){
  
  #iteration number
  iter = iter + 1
  
  for (iz in 1:nz){
    
    for (ik in 1:nk){      

      # Calculate c for each k in each possible state for the economy
      
      c = z[iz] * k[ik]^alpha + ( 1 - delta ) * k[ik] - k
      c = pmax(c, 1e-8)                                     # preventing c from being negative
      
      v = utility(c, gamma) + beta * zprob[iz,] %*% V
      
      newV[iz, ik] = max(v)

      g[iz, ik] = k[ which( v == max(v), arr.ind = TRUE )[2] ]
      
     }
    
  }

  d = norm(newV-V)
  V = newV
  
  print(d)
    
}


#### Graphs ####

# Policy functions

library(ggplot2)
library(latex2exp)

data = data.frame( t(g) , t(k) )
colnames(data) = c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "k")

ggplot(data) + geom_line(aes(x = k, y = p1), size = 0.8) + 
  geom_line(aes(x = k, y = p2), color = "gray", size = 2) + 
  geom_line(aes(x = k, y = p3), color = "green", size = 2) + 
  geom_line(aes(x = k, y = p4), color = "blue", size = 2) + 
  geom_line(aes(x = k, y = p5), color = "red", size = 2) + 
  geom_line(aes(x = k, y = p6), color = "yellow", size = 2) + 
  geom_line(aes(x = k, y = p7), color = "brown", size = 2) + 
      theme_bw() + 
  xlab(TeX("$k_{t}$"))  + ylab(TeX("$k_{t+1}$")) + theme(aspect.ratio=1) + 
  ggtitle("Policy Functions") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=24) ) 


#### Simulation ####

# Simulated variables

zsim = rep(NA, simul + 1)
ksim = rep(NA, simul + 1)
ysim = rep(NA, simul)
csim = rep(NA, simul)
ksim[1] = k_ss

# Simulate the stochastic Markov Process for z

markov <- function(z, P, s0, n){
  
    "
    Simulates the Markov process
    Returns a data frame with the simulated states and the state variable for n periods
    
    Parameters
    ----------
     z: state space of the discretized process
     P: markov transition matrix where P[i, j] is the probability of transitioning from x[i] to x[j]
    y0: inital value of the variable
    s0: initial state
     n: number of points
    "
    y         = rep( NA, n )     # state variable
    y[1]      = z[s0]            # initial value
    states    = rep( NA, n )     # simulated states
    states[1] = s0
    
    for (t in 2:n){
    
    s1          = states[t-1]                        # previous state
    p           = P[s1, ]                            # probability of transition given previous state
    state       = which( rmultinom( 1, 1, p ) == 1)  # new state
    states[ t ] = state  
    y[ t ]      = z[ state ]  
    
    }
    
    return( data.frame( states, z=y ) )
}

s0              = 4                                  # initial state (z = 1)

simulation = markov(z, zprob, s0, simul)

zsim    = simulation$z
zstate  = simulation$states

# Capital accumulation

for (t in 1:simul) {
  
  ksim[t+1] = g[ zstate[ t + 1 ],  which( k == ksim[ t ], arr.ind = TRUE )[2]  ]
  
}

ysim = zsim[1:simul-1] * ksim[1:simul-1]^alpha
csim = ysim + ( 1- delta ) * ksim[2:simul] -  ksim[1:simul-1]


# Output and consumption

data = data.frame( t = seq(1:(simul-1)) , 
                   y = log( ysim / y_ss ) * 100 , 
                   z = log( zsim[1:simul-1]  / z_ss ) * 100 ,
                   k = log( ksim[1:simul-1] / k_ss ) * 100,
                   c = log( csim / c_ss ) * 100)

ggplot(data) + geom_line(aes(x = t, y = z ), size = 0.8, color = "blue" ) + 
  theme_bw() + 
  xlab(TeX(""))  + ylab(TeX("$\\hat{z}_t$")) + theme(aspect.ratio=1) + 
  ggtitle("Technology") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=12) ) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggplot(data) + geom_line(aes(x = t, y = y ), size = 0.8, color = "blue" ) + 
  theme_bw() + 
  xlab(TeX(""))  + ylab(TeX("$\\hat{y}_t$")) + theme(aspect.ratio=1) + 
  ggtitle("Output") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=12) ) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggplot(data) + geom_line(aes(x = t, y = k ), size = 0.8, color = "blue" ) + 
  theme_bw() + 
  xlab(TeX(""))  + ylab(TeX("$\\hat{k}_t$")) + theme(aspect.ratio=1) + 
  ggtitle("Capital") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=12) ) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggplot(data) + geom_line(aes(x = t, y = c ), size = 0.8, color = "blue" ) + 
  theme_bw() + 
  xlab(TeX(""))  + ylab(TeX("$\\hat{c}_t$")) + theme(aspect.ratio=1) + 
  ggtitle("Consumption") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=12) ) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))


# Standard deviation

desv <- matrix(NA, nrow = length(data[2:length(data)]), ncol = 1)
rownames(desv) <- c("Output", "Technology", "Capital", "Consumption" )
colnames(desv) <- c("sd / sd(y)")

for (i in 1:length( desv ) ){
  
  desv[i] = apply( data[2:length(data)][i], 2, sd, na.rm = TRUE)
  desv[i] = desv[i] / apply( data[2:length(data)][1], 2, sd, na.rm = TRUE)
  
}

desv

# Correlation

library(kableExtra)

kable( cor(data[2:length(data)]), format = "latex", booktabs = TRUE)


