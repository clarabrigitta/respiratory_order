# This script makes the interaction epidemic profiles from the PLOS Pathogens paper
# Influenza Interaction with Cocirculating Pathogens, and Its Impact on Surveillance, Pathogenesis and Epidemic Profile: 
#  A Key Role for Mathematical Modeling
# By L Opatowski, M Baguelin, RM Eggo

# Script by RM Eggo. r.eggo@lshtm.ac.uk

# load packages required. If you dont already have them, use install.packages("deSolve")
library(deSolve)
library(shape)

# add the location of the folder where the figure should be made
your.directory <- "~/Desktop/interaction_disruption"

# Make a PDF of the figure
pdf(paste0(your.directory, "Epidemics_Interaction_Figure.pdf"), height=8, width=8, useDingbats = FALSE)

#### BACTERIAL EXAMPLE with enhancement
# define model 
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N = S + I1 + I2 + I12 + I21 + R1 + R2 + R12 + R21 + C12 + C21 # sum all groups
    
    # people getting 1 from single 1s, people getting 1 from co-infected 12s, people getting 1 from coinfected 21s
    # people getting 2 from single 2s, people getting 2 from coinfected 12s, people getting 2 from coinfected 21s
    dS <- (- b1*S*(I1/N + C21/N) - b1*sigma2*S*(I21/N + I12/N) # infections to I1
           - b2*S*(I2/N + C12/N) - b2*sigma1*S*(I12/N + I21/N) ) # infections to I2
    
    # infecteds
    # single infections
    dI1 <- ( + b1*S*(I1/N + C21/N) + b1*sigma2*S*(I21/N + I12/N) # arriving from S
            - delta2*(b2*I1*(I2/N + C12/N) + b2*sigma1*I1*(I12/N + I21/N)) # leaving to double infection
            - gamma1*I1)
    dI2 <- ( + b2*S*(I2/N + C12/N) + b2*sigma1*S*(I12/N + I21/N) # arriving from S
             - delta1*(b1*I2*(I1/N + C21/N) + b1*sigma2*I2*(I21/N + I12/N)) # leaving to double infection
            - gamma2*I2)
    # double infections
    dI12 <- ( + delta2*(b2*I1*(I2/N + C12/N) + b2*sigma1*I1*(I12/N + I21/N))# infections from I1
             - gamma12*I12 )
    dI21 <- ( + delta1*(b1*I2*(I1/N + C21/N) + b1*sigma2*I2*(I21/N + I12/N)) #infections from I2
             - gamma21*I21)
    # consecutively infected
    dC12 <- ( + b2*R1*(I2/N + C12/N) + b2*sigma1*R1*(I12/N + I21/N)  # getting 2 after youve had 1
              - gamma2*C12 )
    dC21 <- ( + b1*R2*(I1/N + C21/N) + b1*sigma2*R2*(I21/N + I12/N) # getting 1 after youve had 2
              - gamma1*C21 )
    
    # recovereds
    # single infections
    dR1 <- gamma1*I1 - b2*R1*(I2/N + C12/N) - b2*sigma1*R1*(I12/N + I21/N) 
    dR2 <- gamma2*I2 - b1*R2*(I1/N + C21/N) - b1*sigma2*R2*(I21/N + I12/N)
    # double infections
    dR12 <- gamma12*I12 + gamma2*C12 
    dR21 <- gamma21*I21 + gamma1*C21
    
    return(list(c(dS, dI1, dI2, dI12, dI21, dC12, dC21, dR1, dR2, dR12, dR21)))   #, dItot1, dItot2
  })
}

# set the constant values
N <- 10000 # total population
times <- seq(0, 100, by = 0.5)   

# initial conditions (NB higher I2 than viral example)
init <- c(S = 9499 ,
          I1 = 1, I2 = 500, I12 = 0, I21 = 0,
          C12 = 0, C21 = 0,
          R1 = 0, R2 = 0, R12 = 0, R21 = 0) 

# scenario 1 parameters
beta1 <- 1 # transmissibility pathogen 1
beta2 <- 0.05 # transmissibility pathogen 2
sigma1 <- 1 # change in infectiousness of path 2 in those who have path 1
sigma2 <- 1 # change in infectiousness of path 1 in those who have path 2
delta1 <- 1 # change in infection rate of path 1 in those who have path 2
delta2 <- 1 # change in infection rate of path 2 in those who have path 1
gamma1 <- 0.25 # recovery from path 1
gamma2 <- 0.05 # recovery from path 2
parameters1 <- list(b1 = beta1, b2 = beta2,
                   sigma1 = sigma1, sigma2 = sigma2,
                   delta1 = delta1, delta2 = delta2,
                   gamma1 = gamma1, gamma2 = gamma2, 
                   gamma12 = gamma2, gamma21 = gamma2) 
# run the ODE for scenario 1
out1 <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters1)) 

# scenario 2 parameters
beta1 <- 1 # transmissibility pathogen 1
beta2 <- 0.05 # transmissibility pathogen 2
sigma1 <- 4 # change in infectiousness of path 2 in those who have path 1
sigma2 <- 1 # change in infectiousness of path 1 in those who have path 2
delta1 <- 1 # change in infection rate of path 1 in those who have path 2
delta2 <- 1 # change in infection rate of path 2 in those who have path 1
gamma1 <- 0.25
gamma2 <- 0.05
parameters2 <- list(b1 = beta1, b2 = beta2,
                   sigma1 = sigma1, sigma2 = sigma2,
                   delta1 = delta1, delta2 = delta2,
                   gamma1 = gamma1, gamma2 = gamma2, 
                   gamma12 = gamma2, gamma21 = gamma2) 
# run the ODE for scenario 2
out2 <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters2)) 

# scenario 3 parameters
beta1 <- 1 # transmissibility pathogen 1
beta2 <- 0.05 # transmissibility pathogen 2
sigma1 <- 4 # change in infectiousness of path 2 in those who have path 1
sigma2 <- 2 # change in infectiousness of path 1 in those who have path 2
delta1 <- 1 # change in infection rate of path 1 in those who have path 2
delta2 <- 1 # change in infection rate of path 2 in those who have path 1
gamma1 <- 0.25
gamma2 <- 0.05
parameters3 <- list(b1 = beta1, b2 = beta2,
                    sigma1 = sigma1, sigma2 = sigma2,
                    delta1 = delta1, delta2 = delta2,
                    gamma1 = gamma1, gamma2 = gamma2, 
                    gamma12 = gamma2, gamma21 = gamma2) 
# run the ODE for scenario 3
out3 <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters3)) 

#### MAKE PLOTS of 3 bacterial scenarios
# vector for colors
cols <- c("black", #S
          "firebrick3","dodgerblue2", #I1, I2
          "mediumpurple1","mediumpurple4", #I12, I21
          "grey45", "grey55", "grey65", "grey75") 

# layout the plot
layout(matrix(c(4, 4, 2, 2,
                1, 1, 2, 2,
                1, 1, 3, 3,
                5, 5, 3, 3,
                9, 9, 7, 7,
                6, 6, 7, 7,
                6, 6, 8, 8,
                10, 10, 8, 8), nrow=8, byrow=TRUE))
#layout.show(n=10)

par( mar=c(2,4,0.5,0.5), las=1, lwd=1.5, mgp=c(3,0.5,0))
xlim <- 40
ylim <- 0.52

# main plot
plot(out1$time, (out1$I1 + out1$I12 + out1$I21 + out1$C21)/N, type = "l", col=cols[2], xlab="time", ylab="", bty = "l", ylim=c(0, ylim), xlim=c(0,xlim), tck=-0.027)  #plot for infectious
lines(out1$time, ((out1$I2 + out1$I12 + out1$I21 + out1$C12)/N), type = "l", col=cols[3], xlab="time", ylab="", bty ="l")  #plot for recovered

# add arrows and legends
x0 <- out1$time[which.max((out1$I1 + out1$I12 + out1$I21 + out1$C21)/N)]
y0 <- max((out1$I1 + out1$I12 + out1$I21 + out1$C21)/N)
# vertical arrow
arrows(x0=x0, y0=0, x1=x0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=x0, y0=0, angle = -90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
# horizontal arrow
arrows(x0=x0, y0=y0, x1=0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=-0.5, y0=y0, angle = 180, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
mtext(side=2, line=2.6, text="Prevalence", las=3, cex=0.9)
mtext(side=1, line=2, text="Time", cex=0.9)
legend("topleft", "No interaction", bty="n")
legend(legend=c("Influenza", "Enhancing pathogen"), col=cols[2:3], lwd=2, bty="n", x=0, y=0.7, xpd=NA)


# second plot
plot(out2$time, (out2$I1 + out2$I12 + out2$I21 + out2$C21)/N, type = "l", col=cols[2], xlab="time", ylab="", bty = "l", ylim=c(0, ylim), xlim=c(0,xlim), tck=-0.027)  #plot for infectious
lines(out2$time, ((out2$I2 + out2$I12 + out2$I21 + out2$C12)/N), type = "l", col=cols[3], xlab="time", ylab="number", bty ="l")  #plot for recovered
y0 <- max((out2$I2 + out2$I12 + out2$I21 + out2$C12)/N)
arrows(x0=10, y0=0.05, x1=45, y1=0.05, lwd=2, lty=3, col=cols[3], length=0, angle=20, code=1)
arrows(x0=40, y0=0.05, x1=40, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=40, y0=y0-0.02, angle = 90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
         arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
Arrowhead(x0=40, y0=0.05+0.02, angle = -90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
legend(x=-2, y=0.57, c("4-fold increase in bacterial transmission during coinfection"), bty="n", xpd=NA)

# third plot
plot(out3$time, (out3$I1 + out3$I12 + out3$I21 + out3$C21)/N, type = "l", col=cols[2], xlab="time", ylab="", bty = "l", ylim=c(0, ylim), xlim=c(0,xlim), tck=-0.027)  #plot for infectious
lines(out3$time, ((out3$I2 + out3$I12 + out3$I21 + out3$C12)/N), type = "l", col=cols[3], xlab="time", ylab="number", bty ="l")  #plot for recovered

x0 <- out3$time[which.max((out3$I1 + out3$I12 + out3$I21 + out3$C21)/N)]
y0 <- max((out3$I1 + out3$I12 + out3$I21 + out3$C21)/N)
# vertical arrow
arrows(x0=x0, y0=0, x1=x0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=x0, y0=0, angle = -90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
# horizontal arrow
arrows(x0=x0, y0=y0, x1=0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=-0.5, y0=y0, angle = 180, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
legend(x=-2, y=0.57, c("4-fold increase in bacterial transmission during coinfection",
                     "plus 2-fold increase in influenza transmission"), bty="n", xpd=NA)

# Add in label
plot.new()
legend("topleft", "A", bty="n", cex=3.5, inset=-0.25, xpd=NA)
plot.new()

##################################################################
### VIRAL EXAMPLE
# The only difference from the bacterial example is the introduction time. 
# the value "test" adds an I2 individual at time = 5. See the end of dI2
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    test = if(floor(time)==5) { test = 1 } else {test = 0}
    
    N = S + I1 + I2 + I12 + I21 + R1 + R2 + R12 + R21 + C12 + C21 # sum all groups
    
    # people getting 1 from single 1s, people getting 1 from co-infected 12s, people getting 1 from coinfected 21s
    # people getting 2 from single 2s, people getting 2 from coinfected 12s, people getting 2 from coinfected 21s
    dS <- (- b1*S*(I1/N + C21/N) - b1*sigma2*S*(I21/N + I12/N) # infections to I1
           - b2*S*(I2/N + C12/N) - b2*sigma1*S*(I12/N + I21/N) ) # infections to I2
    
    # infecteds
    # single infections
    dI1 <- ( + b1*S*(I1/N + C21/N) + b1*sigma2*S*(I21/N + I12/N) # arriving from S
             - delta2*(b2*I1*(I2/N + C12/N) + b2*sigma1*I1*(I12/N + I21/N)) # leaving to double infection
             - gamma1*I1)
    dI2 <- ( + b2*S*(I2/N + C12/N) + b2*sigma1*S*(I12/N + I21/N) # arriving from S
             - delta1*(b1*I2*(I1/N + C21/N) + b1*sigma2*I2*(I21/N + I12/N)) # leaving to double infection
             - gamma2*I2 + test)
    # double infections
    dI12 <- ( + delta2*(b2*I1*(I2/N + C12/N) + b2*sigma1*I1*(I12/N + I21/N))# infections from I1
              - gamma12*I12 )
    dI21 <- ( + delta1*(b1*I2*(I1/N + C21/N) + b1*sigma2*I2*(I21/N + I12/N)) #infections from I2
              - gamma21*I21)
    # consecutively infected
    dC12 <- ( + b2*R1*(I2/N + C12/N) + b2*sigma1*R1*(I12/N + I21/N)  # getting 2 after youve had 1
              - gamma2*C12 )
    dC21 <- ( + b1*R2*(I1/N + C21/N) + b1*sigma2*R2*(I21/N + I12/N) # getting 1 after youve had 2
              - gamma1*C21 )
    
    # recovereds
    # single infections
    dR1 <- gamma1*I1 - b2*R1*(I2/N + C12/N) - b2*sigma1*R1*(I12/N + I21/N) 
    dR2 <- gamma2*I2 - b1*R2*(I1/N + C21/N) - b1*sigma2*R2*(I21/N + I12/N)
    # double infections
    dR12 <- gamma12*I12 + gamma2*C12 
    dR21 <- gamma21*I21 + gamma1*C21
    
    return(list(c(dS, dI1, dI2, dI12, dI21, dC12, dC21, dR1, dR2, dR12, dR21)))   #, dItot1, dItot2
  })
}

# set constants for 3 scenarios
N <- 10000
times <- seq(0, 100, by = 0.5)   

# initial conditions
init <- c(S = 9999 ,
          I1 = 1, I2 = 0, I12 = 0, I21 = 0,
          C12 = 0, C21 = 0,
          R1 = 0, R2 = 0, R12 = 0, R21 = 0) 

# scenario 1 parameters
beta1 <- 1 # transmissibilities
beta2 <- 1
sigma1 <- 1 # change in infectiousness of path 2 in those who have path 1
sigma2 <- 1 # change in infectiousness of path 1 in those who have path 2
delta1 <- 1 # change in infection rate of path 1 in those who have path 2
delta2 <- 1 # change in infection rate of path 2 in those who have path 1
gamma1 <- 0.25
gamma2 <- 0.25
parameters1 <- list(b1 = beta1, b2 = beta2,
                    sigma1 = sigma1, sigma2 = sigma2,
                    delta1 = delta1, delta2 = delta2,
                    gamma1 = gamma1, gamma2 = gamma2, 
                    gamma12 = gamma2, gamma21 = gamma2) 
# run scenario 1 ODE
out1 <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters1)) 

# scenario 2 parameters
beta1 <- 1 # transmissibilities
beta2 <- 1
sigma1 <- 1 # change in infectiousness of path 2 in those who have path 1
sigma2 <- 1 # change in infectiousness of path 1 in those who have path 2
delta1 <- 1 # change in infection rate of path 1 in those who have path 2
delta2 <- 0.6 # change in infection rate of path 2 in those who have path 1
gamma1 <- 0.25
gamma2 <- 0.25
parameters2 <- list(b1 = beta1, b2 = beta2,
                    sigma1 = sigma1, sigma2 = sigma2,
                    delta1 = delta1, delta2 = delta2,
                    gamma1 = gamma1, gamma2 = gamma2, 
                    gamma12 = gamma2, gamma21 = gamma2) 
# run scenario 2 ODE
out2 <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters2)) 

# scenario 3 parameters
beta1 <- 1 # transmissibilities
beta2 <- 1
sigma1 <- 1 # change in infectiousness of path 2 in those who have path 1
sigma2 <- 1 # change in infectiousness of path 1 in those who have path 2
delta1 <- 1 # change in infection rate of path 1 in those who have path 2
delta2 <- 0.1 # change in infection rate of path 2 in those who have path 1
gamma1 <- 0.25
gamma2 <- 0.25
parameters3 <- list(b1 = beta1, b2 = beta2,
                    sigma1 = sigma1, sigma2 = sigma2,
                    delta1 = delta1, delta2 = delta2,
                    gamma1 = gamma1, gamma2 = gamma2, 
                    gamma12 = gamma2, gamma21 = gamma2) 
# run scenario 3 ODE
out3 <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters3)) 

# Plot the 3 viral scenarios
par( mar=c(2,4,0.5,0.5), las=1, lwd=1.5, mgp=c(3,0.5,0))
xlim <- 40
plot(out1$time, (out1$I1 + out1$I12 + out1$I21 + out1$C21)/N, type = "l", col=cols[2], xlab="time", ylab="", bty = "l", ylim=c(0, ylim), xlim=c(0,xlim), tck=-0.027)  #plot for infectious
lines(out1$time[11:100], ((out1$I2 + out1$I12 + out1$I21 + out1$C12)/N)[11:100], type = "l", col=cols[3], xlab="time", ylab="", bty ="l")  #plot for recovered

x0 <- out1$time[which.max((out1$I2 + out1$I12 + out1$I21 + out1$C12)/N)]
y0 <- max((out1$I2 + out1$I12 + out1$I21 + out1$C12)/N)
#vertical arrow
arrows(x0=x0, y0=0, x1=x0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=x0, y0=0, angle = -90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
# horizontal arrow
arrows(x0=x0, y0=y0, x1=0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=-0.5, y0=y0, angle = 180, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
# labels
mtext(side=2, line=2.6, text="Prevalence", las=3, cex=0.9)
arrows(x0 = 5, y0 = -0.05, x1 = 5, y1 = -0.1, length = 0.1, angle = 20,
       code = 1, xpd=NA, col=cols[3])
text(x=5, y=-0.14, "introduction of \ncompeting pathogen", xpd=NA)
mtext(side=1, line=2, text="Time", cex=0.9)
legend("topleft", "No interaction", bty="n")
legend(legend=c("Influenza", "Competing pathogen"), col=cols[2:3], lwd=2, bty="n", x=0, y=0.7, xpd=NA)

# second plot
plot(out2$time, (out2$I1 + out2$I12 + out2$I21 + out2$C21)/N, type = "l", col=cols[2], xlab="time", ylab="", bty = "l", ylim=c(0, ylim), xlim=c(0,xlim), tck=-0.027)  #plot for infectious
lines(out2$time[11:100], ((out2$I2 + out2$I12 + out2$I21 + out2$C12)/N)[11:100], type = "l", col=cols[3], xlab="time", ylab="number", bty ="l")  #plot for recovered

x0 <- out2$time[which.max((out2$I2 + out2$I12 + out2$I21 + out2$C12)/N)]
y0 <- max((out2$I2 + out2$I12 + out2$I21 + out2$C12)/N)
#vertical arrow
arrows(x0=x0, y0=0, x1=x0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=x0, y0=0, angle = -90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
# horizontal arrow
arrows(x0=x0, y0=y0, x1=0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=-0.5, y0=y0, angle = 180, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
legend(x=-2, y=0.54, "Competition decreases chance of coinfection by 50%", bty="n", xpd=NA)

# third plot
plot(out3$time, (out3$I1 + out3$I12 + out3$I21 + out3$C21)/N, type = "l", col=cols[2], xlab="time", ylab="", bty = "l", ylim=c(0, ylim), xlim=c(0,xlim), tck=-0.027)  #plot for infectious
lines(out3$time[11:100], ((out3$I2 + out3$I12 + out3$I21 + out3$C12)/N)[11:100], type = "l", col=cols[3], xlab="time", ylab="number", bty ="l")  #plot for recovered
x0 <- out3$time[which.max((out3$I2 + out3$I12 + out3$I21 + out3$C12)/N)]
y0 <- max((out3$I2 + out3$I12 + out3$I21 + out3$C12)/N)
#vertical arrow
arrows(x0=x0, y0=0, x1=x0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=x0, y0=0, angle = -90, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
# horizontal arrow
arrows(x0=x0, y0=y0, x1=0, y1=y0, lwd=2, lty=3, col=cols[9], length=0, angle=20, code=1)
Arrowhead(x0=-0.5, y0=y0, angle = 180, arr.length = 0.2, arr.width = 0.15, arr.adj = 0.5, 
          arr.type = "triangle", lcol = cols[9], lty = 1, arr.lwd = 1.5, npoint = 5)
legend(x=-2, y=0.54, "Competition decreases chance of coinfection by 90%", bty="n", xpd=NA)

# add label
plot.new()
legend("topleft", "B", bty="n", cex=3.5, inset=-0.25, xpd=NA)
plot.new()

# close plot
dev.off()

##### END #######################################################
