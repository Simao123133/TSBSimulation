n = 2

info = 0
iflag = 0
nfev = 0

iflag = 1
call fcn(n,x,fvec,iflag)
nfev = 1
fnorm = enorm(n,fvec)

#determine the number of calls to fcn needed to compute the jacobian matrix.
msum = min0(ml+mu+1,n)

#initialize iteration counter and monitors.

iter = 1
ncsuc = 0
ncfail = 0
nslow1 = 0
nslow2 = 0

#beginning of the outer loop.
iflag = 2
call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
nfev = nfev + msum

#compute the qr factorization of the jacobian.
call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)

#on the first iteration and if mode is 1, scale according to the norms of the columns of the initial jacobian.
