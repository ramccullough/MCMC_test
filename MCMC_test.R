alpha_true=2
beta_true=1
#sd_true=1
n_samples=48
n_params=2

acc_lower=0.35
acc_upper=0.45

#create sample data by adding normally distributed noise
x=seq(from=0, to=12, by=0.25)
y=alpha_true*x+beta_true+rnorm(length(x), 0, 1)


N=100000 #number of steps in chain
#alpha_values<-matrix(0, N, 1)
#beta_values<-matrix(0, N, 1)
mchain=matrix(0, N, n_params) #set up MC matrix for parameters 

alpha_init<-rnorm(1, 2, 1) #initial guess for alpha (slope)
beta_init<-rnorm(1, 1, 1)
#sigma_init<-rnorm(1, 1, 1)

#alpha_values[, 1]<-alpha_init
#beta_values[, 2]<-beta_init

mchain[1, 1]<-alpha_init
mchain[1, 2]<-beta_init
#mchain[1, 3]<-sigma_init

#define prior distribution
prior=function(param){
   dunif(param, min=-1000, max=1000, log=TRUE)
}

#define proposal distribution to draw guesses from
proposal=function(param, s){
  return(rnorm(1, mean=param, sd=s))
}

#initialise sd
s=matrix(1, 1, n_params)

#initialise accepted values count
count=matrix(0, 1, n_params)

for (t in 2:N){
  for (i in 1:n_params) {
    #tuning sd for proposal distribution
    
    if (t<=N/2 && t%%100==0){
      acc_rate=count[, i]/t
      if (acc_rate<acc_lower){
        scale=acc_lower/acc_rate
        s[, i]=s[, i]/scale}
      if (acc_rate>acc_upper){
        scale=acc_rate/acc_upper
        s[, i]=s[, i]*scale}}
    
    #choose value from proposal density   
    k=proposal(mchain[t-1, i], s[, i])
    
    #alpha likelihood 
    if (i==1){Likelihood_old=sum(dnorm(y, mean=mchain[t-1, i]*x+mchain[t-1, i+1], sd=1, log=T))
    Likelihood_new=sum(dnorm(y, mean=k*x+mchain[t-1, i+1], sd=1, log=T))}
    
    #beta likelihood
    if (i==2){Likelihood_old=sum(dnorm(y, mean=mchain[t, i-1]*x+mchain[t-1, i], sd=1, log=T)) 
    Likelihood_new=sum(dnorm(y, mean=mchain[t, i-1]*x+k, sd=1, log=T))}
    
    #sigma likelihood
    # if (i==3){Likelihood_old=sum(dnorm(y, mean=mchain[t, i-2]*x+mchain[t, i-1], sd=log(mchain[t-1, i]), log=T))
    #Likelihood_new=sum(dnorm(y, mean=mchain[t, i-2]*x+mchain[t, i-1], sd=log(k), log=T))}
    
    #calculate priors
    Prior_new=prior(k)
    Prior_old=prior(mchain[t-1, i])
    
    #calculate posteriors
    Posterior_new=Likelihood_new+Prior_new
    Posterior_old=Likelihood_old+Prior_old
    
    #calculate acceptance probability
    #a=min(1, exp(Posterior_new-Posterior_old))
    log_a=min(0, Posterior_new-Posterior_old)
    
    #choose random value from uniform distribution
    u=runif(1, 0, 1)
    log_u=log(u)
    
    if (log_a>log_u) {
      mchain[t, i]=k
      count[, i]=count[, i]+1
      } 
    else {
      mchain[t, i]=mchain[t-1, i]
    }
    #if (log_a>log_u) {
    #if (i==1 || i==2){
    
     # mchain[t, i]=k
    #  count[, i]=count[, i]+1}
    #if (i==3){
    #mchain[t, i]=exp(k)
    #}
    #} 
    #else {
      #mchain[t, i]=mchain[t-1, i]
    #} 
  }
}

#trace plots
plot(1:N, mchain[, 1], type="l", xlab="Iteration", ylab="Alpha")
plot(1:N, mchain[, 2], type="l", xlab="Iteration", ylab="Beta")

#remove burn-in iterations
data=mchain[50001:N, ]
alpha_post=data[, 1]
beta_post=data[, 2]

#trace plots without burn-in
plot(50001:N, alpha_post, type="l", xlab="Iteration", ylab="Alpha")
plot(50001:N, beta_post, type="l", xlab="Iteration", ylab="Beta")

hist(alpha_post, freq=F, xlab="alpha")
curve(dnorm(x,mean=2, sd=s[, 1]), add=TRUE, lwd=2)

hist(beta_post, freq=F, xlab="beta")
curve(dnorm(x,mean=1, sd=s[, 2]), add=TRUE, lwd=2)

#acceptance rates
acceptance_alpha=count[, 1]/N
acceptance_beta=count[, 2]/N

