``` r
n = 200    # taille de la population
T = 1000   # nombre de pas de temps
V = matrix(0,n,T)    # V = charge virale
S = matrix(1,n,T)    # matrice des S
I = matrix(0,n,T)    # matrice des I
S[1,1:T]=0  # la personne 1 est infectée
I[1,1]=1
R = matrix(0,n,T)    # matrice des R
pop = 1:n
lambda = .1    # nombre de contacts, par date (Poisson)

set.seed(1)    # simulation d'une courbe de charge virale
               # je prends une loi Gamma
               # on est contagieux seulement si le charge virale dépasse 1
charge = dgamma((1:T)/100,3,20)*runif(1)
V[1,1:T] = ifelse(charge>.1,charge,0)
I[1,1:T] = ifelse(charge>.1,1,0)
R[1,1:T] = 1-S[1,1:T]-I[1,1:T]

# on itère, date par date
for(t in 2:T){
  i = which(I[,t] == 1) # on prend toutes les personnes infectées
  for(idx in i){
    nb=rpois(1,lambda)  # on tire le nombre de contacts que la personne aura à la date t (dans toute la population)
    if(nb>0){
    indiv_contact = sample(pop,size=nb)
    proba = pmin(V[idx,t],1)
    contamine = rbinom(nb,1,proba) # elle contamine si la charge virale est assez`grande, et si le contact est un S
    indiv_contamine = indiv_contact[contamine&S[indiv_contact,t]]
   
    if(length(indiv_contamine)>0){ # pour chaque contact contaminé, on tire sa charge virale, etc
    for(j in indiv_contamine)
      charge = dgamma((1:T)/100,3,20)[1:(T-t+1)]*runif(1)
      signe = c(1,diff(charge))>0
      V[j,t:T] = ifelse(charge>.1,charge,0)
      S[j,t:T] = ifelse((charge>.1)|((signe==0)&max(charge)>.1),0,1)
      I[j,t:T] = ifelse(charge>.1,1,0)
      R[j,t:T] = 1-S[j,t:T]-I[j,t:T]
    }}
  }}
plot(1:T,apply(S,2,mean),type="l",col="blue",ylim=0:1)
lines(1:T,apply(I,2,mean),type="l",col="red")
```
