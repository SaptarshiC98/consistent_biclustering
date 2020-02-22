############__________ l_2 distance_________________###############
S_new=function(x,y){
  a=abs(x)
  if(a>y){
    return(sign(x)*(a-y))
  }
  else{
    return(0)
  }
}
sparse_bc=function(X,K,R,lambda,tmax=100){
  Cs <- kmeans(X, K,nstart = 20)$cluster
  Dr <- kmeans(t(X),R,nstart = 20)$cluster
  MU=matrix(0,K,R)
  dk=numeric(K)
  dr=numeric(R)
  n=dim(X)[1]
  p=dim(X)[2]
  m=X
  # UPDATE MU
  for(t in 1:tmax){
    for(k in 1:K){
      I=which(Cs==k)
      for(r in 1:R){
        J=which(Dr==r)
        MU[k,r]=S_new(sum(X[I,J]),lambda)/(length(I)*length(J))
        #cat(c(k,r))
        #cat('\n')
      }
    }
    
    for(i in 1:n){
      for(k in 1:K){
        dk[k]=0
        for(j in 1:p){
          dk[k]=dk[k]+(X[i,j]-MU[k,Dr[j]])^2
        }
      }
      Cs[i]=which.min(dk)
    }
    
    for(j in 1:p){
      for(r in 1:R){
        dr[r]=0
        for(i in 1:n){
          dr[r]=dr[r]+(X[i,j]-MU[Cs[i],r])^2
        }
        Dr[j]=which.min(dr)
      }
    }
  }
  for(i in 1:n){
    for(j in 1:p){
      m[i,j]=MU[Cs[i],Dr[j]]
    }
  }
  return(list(Cs,Dr,MU,m))
}
#############______________ l_1 distance___________###############
Update_mu_l1=function(x,lambda){
  n=length(x)
  d=function(x,y){
    return(sum(abs(x-y)))
  }
  f=function(mu){
    s=0
    for(i in 1:n){
      s=s+d(x[i],mu)
    }
    s=s+lambda*sum(abs(mu))
    return(s)
  }
  a=min(x)
  b=max(x)
  if(a>0){
    a=-1
  }else if(b<0){
    b=1
  }
  I=c(a,b)
  o=optimize(f,I)
  return(o$minimum)
}
li_dist=function(x,y){
  return(sum(abs(x-y)))
}
sparse_bc_l1=function(X,K,R,lambda,tmax=100){
  Cs=truthCs
  Dr=truthDs
  MU=matrix(0,K,R)
  dk=numeric(K)
  dr=numeric(R)
  n=dim(X)[1]
  p=dim(X)[2]
  m=X
  # UPDATE MU
  for(t in 1:tmax){
    for(k in 1:K){
      I=which(Cs==k)
      for(r in 1:R){
        J=which(Dr==r)
        xx=c(X[I,J])
        MU[k,r]=Update_mu_l1(xx,lambda)
        #cat(c(k,r))
        #cat('\n')
      }
    }
    
    for(i in 1:n){
      for(k in 1:K){
        dk[k]=0
        for(j in 1:p){
          dk[k]=dk[k]+abs(X[i,j]-MU[k,Dr[j]])
        }
      }
      Cs[i]=which.min(dk)
    }
    
    for(j in 1:p){
      for(r in 1:R){
        dr[r]=0
        for(i in 1:n){
          dr[r]=dr[r]+abs(X[i,j]-MU[Cs[i],r])
        }
        Dr[j]=which.min(dr)
      }
    }
    if(t%%10==0){
      cat(t)
      cat('\n')
    }
  }
  for(i in 1:n){
    for(j in 1:p){
      m[i,j]=MU[Cs[i],Dr[j]]
    }
  }
  return(list(Cs,Dr,MU,m))
}
##########_______Huber Loss______############

d_huber=function(x,y,delta=1){
  a=x-y
  if(abs(a)<delta){
    return(a^2/2)
  }else{
    return(delta*(abs(a)-delta/2))
  }
}

Update_mu_huber=function(x,lambda){
  n=length(x)
  f=function(mu){
    s=0
    for(i in 1:n){
      s=s+d_huber(x[i],mu)
    }
    s=s+lambda*abs(mu)
    return(s)
  }
  I=c(min(x),max(x))
  o=optimize(f,I)
  return(o$minimum)
}

sparse_bc_huber=function(X,K,R,lambda,tmax=100){
  Cs <- kmeans(X, K,nstart = 20)$cluster
  Dr <- kmeans(t(X),R,nstart = 20)$cluster
  MU=matrix(0,K,R)
  dk=numeric(K)
  dr=numeric(R)
  n=dim(X)[1]
  p=dim(X)[2]
  m=X
  # UPDATE MU
  for(t in 1:tmax){
    for(k in 1:K){
      I=which(Cs==k)
      for(r in 1:R){
        J=which(Dr==r)
        xx=c(X[I,J])
        MU[k,r]=Update_mu_huber(xx,lambda)
        #cat(c(k,r))
        #cat('\n')
      }
    }
    
    for(i in 1:n){
      for(k in 1:K){
        dk[k]=0
        for(j in 1:p){
          dk[k]=dk[k]+d_huber(X[i,j],MU[k,Dr[j]])
        }
      }
      Cs[i]=which.min(dk)
    }
    
    for(j in 1:p){
      for(r in 1:R){
        dr[r]=0
        for(i in 1:n){
          dr[r]=dr[r]+d_huber(X[i,j],MU[Cs[i],r])
        }
        Dr[j]=which.min(dr)
      }
    }
    if(t%%10==0){
      cat(t)
      cat('\n')
    }
  }
  for(i in 1:n){
    for(j in 1:p){
      m[i,j]=MU[Cs[i],Dr[j]]
    }
  }
  return(list(Cs,Dr,MU,m))
}

#########_______Minkwoski______###########


d_mwk=function(x,y,beta=4){
  return(abs(x-y)^beta)
}

Update_mu_mwk=function(x,lambda,beta=4){
  n=length(x)
  f=function(mu){
    s=0
    for(i in 1:n){
      s=s+d_mwk(x[i],mu,beta)
    }
    s=s+lambda*abs(mu)
    return(s)
  }
  I=c(min(x),max(x))
  o=optimize(f,I)
  return(o$minimum)
}

sparse_bc_mwk=function(X,K,R,lambda,beta=4,tmax=100){
  Cs <- kmeans(X, K,nstart = 20)$cluster
  Dr <- kmeans(t(X),R,nstart = 20)$cluster
  MU=matrix(0,K,R)
  dk=numeric(K)
  dr=numeric(R)
  n=dim(X)[1]
  p=dim(X)[2]
  m=X
  # UPDATE MU
  for(t in 1:tmax){
    for(k in 1:K){
      I=which(Cs==k)
      for(r in 1:R){
        J=which(Dr==r)
        xx=c(X[I,J])
        MU[k,r]=Update_mu_mwk(xx,lambda,beta)
        #cat(c(k,r))
        #cat('\n')
      }
    }
    
    for(i in 1:n){
      for(k in 1:K){
        dk[k]=0
        for(j in 1:p){
          dk[k]=dk[k]+d_mwk(X[i,j],MU[k,Dr[j]])
        }
      }
      Cs[i]=which.min(dk)
    }
    
    for(j in 1:p){
      for(r in 1:R){
        dr[r]=0
        for(i in 1:n){
          dr[r]=dr[r]+d_mwk(X[i,j],MU[Cs[i],r])
        }
        Dr[j]=which.min(dr)
      }
    }
    if(t%%10==0){
      cat(t)
      cat('\n')
    }
  }
  for(i in 1:n){
    for(j in 1:p){
      m[i,j]=MU[Cs[i],Dr[j]]
    }
  }
  return(list(Cs,Dr,MU,m))
}
