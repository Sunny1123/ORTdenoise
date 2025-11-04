library(foreach)
library(doParallel)

getnimg = function(img,sigma,seed=Sys.time())
{
  return(img+matrix(rnorm(nrow(img)*ncol(img),0,sigma),nrow=nrow(img),ncol=ncol(img)))
}

genroot = function(img)
{
  y=c(img)
  x1=list()
  for(i in 1:nrow(img)){
    x1[[i]]=c(1:ncol(img))
  }
  x1=unlist(x1)
  x2=list()
  for(i in 1:ncol(img)){
    x2[[i]]=rep(i,nrow(img))
  }
  x2=unlist(x2)
  X=cbind(rep(1,length(x1)),x1,x2)
  return(list(intensity=y,coord=X,leaf=0))
}

getresval = function(node,w)
{
  x=node$coord
  y=node$intensity
  w = w/sum(w^2)^0.5
  xmod = x%*%w
  fit1=y[xmod>0]
  fit2=y[xmod<=0]
  #res= sum(fit1)^2/(length(fit1)+0.00000000001)+sum(fit2)^2/(length(fit2)+0.0000000000001)
  n=length(y)
  n1=length(fit1)
  n2=length(fit2)
  if(n1<=1 || n2<=1 || n<=1) return(0)
  res=var(y)-(n1/n)*var(fit1)-(n2/n)*var(fit2)
  if(is.na(res)) cat(n1,n2,n)
  return(res)
}

optim_func_gen = function(node)
{
  f = function(w)
  {
    res = getresval(node,w)
    return(-res)
  }
  return(f)
}

getchild = function(node,minsize,minlossimprove=0.01)
{
  if(node$leaf==1)
  {
    return(list(node))
  }
  y=node$intensity
  X=node$coord
  optim_func= optim_func_gen(node)
  w_init=c(-ncol(X[,-1])*mean(X[,-1]),rep(1,ncol(X)-1))
  w_opti=optim(w_init,optim_func,control=list(maxit=1000000))
  w_conv=w_opti$conv
  w_opti=w_opti$par
  # if(w_conv)cat("conv--",w_conv,'\n')
  xmod=X%*%w_opti
  fit1=xmod>0
  fit2=xmod<=0
  lsimp=-optim_func(w_opti)
  child_left= list(intensity = y[fit1],coord = X[fit1,],leaf=0)
  child_right=list(intensity = y[fit2],coord = X[fit2,],leaf=0)
  if((lsimp > minlossimprove) &&
     (length(child_left$intensity)>minsize &&
      length(child_right$intensity)>minsize)
  )
  {
    return(list(child_left,child_right))
  }
  else
  {
    res=node
    res$leaf=1
    return(list(res))
  }
}


getnextgen=function(leaf_nodes,minsize,minlossimprove=0.01)
{
  res=list()
  it=1
  for(node in leaf_nodes)
  {
    if(node$leaf==1)
    {
      res=append(res,list(node))
    }
    else
    {
      temp=getchild(node,minsize,minlossimprove)
      res=append(res,temp)
    }
    it=it+1
  }
  return(res)
}

gettree = function(root,minlossimprove=0.01,minsize=1)
{
  res=list(root)
  it=0
  while(1)
  {
    res1=getnextgen(res,minsize,minlossimprove)
    if(length(res)==length(res1))
    {
      return(res1)
    }
    res=res1
    it=it+1

  }

}




Marknodes = function(tree,img)
{
  res = array(img,dim=c(dim(img),2))
  res[,,1]=img
  for(i in 1:length(tree))
  {
    for(j in 1:nrow(tree[[i]]$coord))
    {
      v=tree[[i]]$coord[j,2:3]
      res[v[1],v[2],2]=i;
    }
  }
  return(res)
}

Rleaftomat = function(tree,img,smband)
{
  # Set up parallel computing
  cl <- makeCluster(detectCores(logical = TRUE)-2)
  registerDoParallel(cl)
  # `%dopar%` <- foreach::`%dopar%`
  # `%do%` <- foreach::`%do%`
  # `%:%` <- foreach::`%:%`
  res = img
  nr=nrow(img)
  nc=ncol(img)
  img = Marknodes(tree,img)
  res = foreach(i= 1:nrow(img[,,1]),.combine=rbind,.inorder=FALSE) %:%
    foreach(j= 1:ncol(img[,,1]),.combine=cbind,.inorder=FALSE) %dopar%
    {
      temp=0
      count=0
      for(p in (-smband):smband)
      {
        loc1=nr-abs(nr-1-abs(i+p-1))
        for(q in (-smband):smband)
        {
          loc2=nc-abs(nc-1-abs(j+q-1))
          if(i==loc1 && j==loc2)
          {
            wt=1
          }
          else {
            wt = exp(-(img[loc1,loc2,1]-img[i,j,1])/sum((c(i,j)-c(loc1,loc2))^2)^0.5)
          }
          temp= temp+img[loc1,loc2,1]*(img[loc1,loc2,2]==img[i,j,2])*wt
          count =  count + (img[loc1,loc2,2]==img[i,j,2])*wt
        }
      }
      return(temp/count)
    }
  return(res)
}

Batchrun = function(n,sigma,tol,img,smband=c(1,2,3,4,5))
{
  # Set up parallel computing
  cl <- makeCluster(min(detectCores(logical = TRUE)-2,n))
  registerDoParallel(cl)
  x1=list()
  for(i in 1:nrow(img)){
    x1[[i]]=c(1:nrow(img))
  }
  x1=unlist(x1)
  x2=list()
  for(i in 1:nrow(img)){
    x2[[i]]=rep(i,nrow(img))
  }
  x2=unlist(x2)
  X=cbind(rep(1,length(x1)),x1,x2)
  #.export=c("gettree","Rleaftomat","getnextgen","getchild","Marknodes","optim_func_gen","getresval","ConvertLeaftoMat","randomForest","JPLLK_surface")
  res=foreach(i = 1:n,.inorder = FALSE,.packages = c("dentree","DRIP","randomForest")) %dopar%
    {
      set.seed(i^2+i+23)
      nimg = img+matrix(rnorm(nrow(img)*ncol(img),0,sigma),nrow=nrow(img),ncol=ncol(img)) # Noisy image
      y=c(nimg)
      root=list(intensity=y,coord=X,leaf=0)
      tree= gettree(root,tol)
      res1 =array(0,dim=c(2,length(smband)))
      row.names(res1)= c("ORT","JPLLK")
      # X=X[,-1]
      # fit=randomForest(X,y,ntree = 25)
      # res1=mean((matrix(predict(fit,X,type="response"),nrow(nimg),ncol(nimg))-img)^2)^0.5
      for(j in 1:length(smband))
      {
        # res1[[j]]=list(ORT=list(),JPLLK=list(),RF=list())
        res1[1,j] = mean((ConvertLeaftoMat(tree,nimg,smband[j])-img)^2)^0.5
        res1[2,j] = mean((JPLLK_surface(nimg,(smband[j]+1),plot = FALSE)$fitted-img)^2)^0.5
        # X=X[,-1]
        # fit=randomForest(X,y,ntree = 25)
        # res1[1,j]=mean((matrix(predict(fit,X,type="response"),nrow(nimg),ncol(nimg))-img)^2)^0.5
      }
      #return(res1)
      return(res1)
    }
  return(res)
}

SSimg = function(img)
{
  resimg = array(0,dim=2*dim(img)-1)
  resimg[2*(1:nrow(img))-1,2*(1:ncol(img))-1] = img
  for(i in 2*(2:nrow(img)-1))
  {
    resimg[i,]=(resimg[i-1,]+resimg[i+1,])/2
  }
  for(j in 2*(2:ncol(img)-1))
  {
    resimg[,j]=(resimg[,j-1]+resimg[,j+1])/2
  }
  return(resimg)
}

edgeMSE=function(img,res,plot=FALSE)
{
  edge=stepEdge(image = img,bandwidth = 4,degree = 0,thresh = 0.2)
  edge1=stepEdge(image = res,bandwidth = 4,degree = 0,thresh = 0.2)
  if(plot)
  {
    image(edge,col=gray(0:256/256))
    image(edge1,col=gray(0:256/256))
  }
  return(dKQ(edge,edge1))
}
batchedge=function(img,sd=c(0.1,0.2,0.3),bw=c(1:5))
{
  cl <- makeCluster(min(detectCores(logical = TRUE),length(bw)))
  registerDoParallel(cl)
  x1=list()
  for(i in 1:nrow(img)){
    x1[[i]]=c(1:nrow(img))
  }
  x1=unlist(x1)
  x2=list()
  for(i in 1:nrow(img)){
    x2[[i]]=rep(i,nrow(img))
  }
  x2=unlist(x2)
  X=cbind(rep(1,length(x1)),x1,x2)
  res=list()
  res=foreach(i = 1:length(bw),.packages = c("dentree","DRIP"))
  {
    res1=matrix(0,nrow=2,ncol=length(sd))
    for(j in 1:length(sd))
    {
      nimg = img+matrix(rnorm(nrow(img)*ncol(img),0,sd[j]),nrow=nrow(img),ncol=ncol(img)) # Noisy image
      y=c(nimg)
      root=list(intensity=y,coord=X,leaf=0)
      tree= gettree(root,tol)
      imgres= ConvertLeaftoMat(tree,nimg,bw[i])
      res1[1,j]=edgeMSE(img,imgres)
      imgJPLLK=JPLLK_surface(img,bw[i]+1)$fitted
      res1[2,j]=edgeMSE(img,imgJPLLK)
    }
    # res[[i]]=res1
    return(res1)
  }
  return(res)
}
