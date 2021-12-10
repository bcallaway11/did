


## distance between logicals
dL <- expand.grid(c(TRUE,FALSE),c(TRUE,FALSE))
x = data.frame(x=dL[,1])
y = data.frame(x=dL[,2])
expect_equal(c(0,1,1,NaN), gower_dist(x = x,y = y))


## distance between factor variables
bands <- c("Grand Magus","Skull Fist")
dF <- expand.grid(bands,bands)
expect_equal(gower_dist(data.frame(x=dF[,1]),data.frame(x=dF[,2])),c(0,1,1,0))

## distance between numerical variables 
dN <- data.frame(x = as.numeric(1:4),y=as.numeric(c(1,1,2,3)))
expect_equal(gower_dist(data.frame(x=dN[,1]),data.frame(x=dN[,2])),c(0,1/3,1/3,1/3))

## distance between character variables
dC <- data.frame(x=letters[1:3],y=letters[3:1],stringsAsFactors=FALSE)
expect_equal(gower_dist( 
  data.frame(x=dC[,1],stringsAsFactors=FALSE)
  , data.frame(x=dC[,2],stringsAsFactors=FALSE)),c(1,0,1))



## multivariate dataset
bands <- c("Grand Magus","Skull Fist")
dL <- expand.grid(c(TRUE,FALSE),c(TRUE,FALSE))
dN <- data.frame(x = as.numeric(1:4),y=as.numeric(c(1,1,2,3)))
dF <- expand.grid(bands,bands)
dM1 <- data.frame(x=dL[,1],y=dF[,1],z=dN[,1])  
dM2 <- data.frame(x=dL[,2],y=dF[,2],z=dN[,2])
expect_equal(gower_dist(x=dM1,y=dM2), c(0,7/9,7/9,1/6))
# check symmetry
expect_equal(gower_dist(dM1,dM2),gower_dist(dM2,dM1))
# not counting NA's in the denominator
dM1[array(c(2,3,4,1,2,3),dim=c(3,2))] <- NA
expect_equal(gower_dist(dM1,dM2), c(0,3/4,3/4,0))
#auto-matching columns
expect_equivalent(gower_dist(women, women[1]),rep(0,nrow(women)))
 


## ignoring column name cases",{
dat1 <- iris[1:6,]
dat2 <- iris[1:6,]
names(dat2) <- tolower(names(dat2))
expect_equal(gower_dist(dat1, dat2, ignore_case=TRUE),rep(0,6))


## recycling
expect_equal(length(gower_dist(x=iris[1,],y=iris)), nrow(iris))
expect_equal(length(gower_dist(x=iris,y=iris[1,])), nrow(iris))
expect_equal(length(gower_dist(x=iris[1:3,],y=iris)), nrow(iris))
expect_equal(length(gower_dist(x=iris,y=iris[1:3,])), nrow(iris))



## weights
expect_error(gower_topn(women, women, weights=-(1:4)))
expect_error(gower_dist(women, women, weights=c(NA,1:3)))

d1 <- women[1,]
d2 <- women[2,]
w <- c(1,2)
r <- sapply(women, function(x) abs(diff(range(x))))

d12 <- (w[1]*abs(d1[1,1]-d2[1,1])/r[1] + w[2]*abs(d1[1,2]-d2[1,2])/r[2])/sum(w)
wom2 <- women
wom2[1:2,] <- wom2[2:1,]
expect_equivalent(gower_dist(women,wom2, weights=w)[1], d12)


## edge cases and exceptions
expect_warning(gower_dist(
  x = data.frame(x=c(1.2,1.2,1.2))
  , y = data.frame(x=c(1.2,1.2,1.2))
))
expect_warning(gower_dist(
  x = data.frame(x=c(1.2,1.2,1.2))
  , y = data.frame(x=c(1.2,1.2,1.3))
  , eps=0.2
))

expect_warning(gower_dist(data.frame(x=rep(1,100)), data.frame(x=1,100)))


expect_error(gower_dist(
  data.frame(a = letters[1:3], stringsAsFactors = TRUE),
  data.frame(a = letters[2:4], stringsAsFactors = TRUE)
))
expect_error(gower_dist(
  data.frame(a = letters[1:3], stringsAsFactors = FALSE),
  data.frame(a = letters[2:4], stringsAsFactors = TRUE)
))

suppressMessages(out <- gower_dist(data.frame(x=1:3), data.frame(y=1:3) ))
expect_true(identical(out,numeric(0)))

suppressMessages(out <- gower_topn(data.frame(x=1:3),data.frame(y=1:3)) )
expect_true( identical(out$distance, matrix(0)[0,0]) )
expect_true( identical(out$index, matrix(0L)[0,0]) )





## Top-n

## gower_topn
d1 <- iris[1:3,]
d2 <- iris[1:7,]
L <- gower_topn(d1,d2,n=4)
expect_equal(length(L),2)
expect_equal(dim(L[[1]]),c(4,3))
expect_equal(dim(L[[1]]),dim(L[[2]]))
expect_equal(L[[1]][1,],1:3)
expect_equal(L[[2]][1,],rep(0,3))

# case where n exceeds nr of records in the lookup table.
L <- gower_topn(d1,d2,n=8)
expect_equal(L[[1]][8,],rep(0,3))
expect_equal(L[[2]][8,],rep(Inf,3))
  



## just to get code-coverage right
dat1 <- data.frame(
  x = as.factor(sample(letters,2000,replace=TRUE))
  ,y = sample(LETTERS,2000,replace=TRUE)
  ,z = as.integer(1:2000)
  ,w = sample(c(TRUE,FALSE),2000,replace=TRUE)
  , stringsAsFactors=FALSE
)
i <- sample(2000)
dat2 <- dat1[i,]
gower_dist(dat1,dat2)



## warning on recycling
expect_warning(gower_dist(iris[1:3,],iris[1:2,]))
expect_warning(gower_dist(iris[1:2,],iris[1:3,]))


## regression tests

## topn w/n=1"
dat <- data.frame(x = c(NA,2,4,5), y = c(6,7,NA,10))
L <- gower_topn(dat[c(1,3),],dat[c(2,4),],n=1)
expect_equivalent(as.vector(L$index),c(1,2))







