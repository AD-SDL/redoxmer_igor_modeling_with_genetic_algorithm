// x is a parameter matrix m x n  ::  n sets of m parameters
// y is a target matrix k x n :: n sets of k targets to fit
// beta minimizes || y - beta * x || , it is m x n matrix
// xx is xx' it is m x m matrix
// yx is yx' it is k x m matrix
// yp = beta * x - predicted values
// res = y-yp - residues
// x and y are first shifted by the means and normalized by variances

//Macro SimpleTestLinReg()
Variable /G yshow=0,r2show=nan
Make /N=0 /O /R yp1,yv1
Make /O line45={-5,-4,-2,-1,0,1,2,3,4,5}
Make /N=1 /O M_product,M_A,M_B
//
Make /N=(3,50) /O /R xv  // 5 points, 2 parameters
Make /N=(2,50) /O /R yv  // 5 points, 1 target
xv = sqrt(p*p+q) 
yv[0][] = 5 * xv[0][q] +   2 * xv[1][q]  + 3 * xv[2][q] +  enoise(8) 
yv[1][] = 8* xv[0][q] -   3 * xv[1][q]  + .2 * xv[2][q] + enoise(4) 
//

ShiftXY(xv,yv)
print LinReg_(xv,yv)
SetYShowProc("",0,"","")
end



function ShiftXY(xv,yv)
Wave xv,yv
Make /N=0 /R /O ax,sx,ay,sy
ShiftVectors(xv,ax,sx,1)
ShiftVectors(yv,ay,sy,1)
end

function ShiftVectors(v,av,sv,opt) // standard form
Wave v,av,sv
Variable opt
Variable m=DimSize(v,0)
Variable n=DimSize(v,1)
Variable i,j,k
Variable z,zx,zxx
for(j=0;j<m;j+=1)
	zx=0; zxx=0; k=0
	for(i=0;i<n;i+=1)
		z=v[j][i]
		if(!numtype(z))
			zx+=z
			zxx+=z*z
			k+=1
		endif
	endfor
	if(k)
		zx/=k
		zxx=zxx/k-zx*zx
		zxx=sqrt(max(0,zxx))
	endif
	av[j]=zx
	sv[j]=zxx
endfor
if(opt)
	v=(v-av[p])/sv[p]
endif
end

function RestoreVectors(v,av,sv)
Wave v,av,sv
	v=v*sv[p]+av[p]
end

function LinRegB_(xvb,yvb,ypb)
Wave xvb,yvb,ypb
Wave beta
Wave M_product
CheckMatrixMultiply(beta,xvb,"LinRegB_")
MatrixMultiply beta,xvb
duplicate /O M_product ypb
Variable /G rmswb=crms(yvb,ypb)
Variable /G r2wb=CalcWilks(yvb,ypb,0) 
return rmswb
end

function LinReg_(xv,yv,yp)
Wave xv,yv,yp
Variable /G lambda
Wave M_product,M_B
MatrixMultiply xv, xv /T
duplicate /O M_product xx // hat matrix
MatrixEigenV /SYM /O xx
xx+=lambda*(p==q) // Ridge parameter added
MatrixMultiply xv, yv /T
duplicate /O M_product yx
MatrixLLS /M=1 xx yx
duplicate /O M_B beta
MatrixTranspose beta
MatrixMultiply beta,xv
duplicate /O M_product yp
Variable /G rmsw=crms(yv,yp)
Variable /G r2w=CalcWilks(yv,yp,1) 
//return r2w // Wilks
return rmsw
end

function Lin_(xv,yv,beta)
Wave xv,yv,beta
Wave M_product,M_B
MatrixMultiply beta,xv
duplicate /O M_product yp
Variable /G rmsw=crms(yv,yp)
Variable /G r2w=CalcWilks(yv,yp,1) 
//return r2w // Wilks
return rmsw
end

function CheckDim(w,point)
Wave w
Variable point
Variable a=DimSize(w,0)
Variable b=DimSize(w,1)
if((!a) || (!b))
	print "zero dimension at point",point
	Abort
endif
end


function LinRegFull(xv,yv,yp) // linear rgression with checks of NaN
Wave xv,yv,yp
Variable m=DimSize(xv,0)
Variable n=DimSize(xv,1)
Variable ny=DimSize(yv,0)
Make /N=(m,m) /O /D xx
Make /N=(m,ny) /O /D yx
Variable i,j,k,z
// X*X' matrix
for(i=0;i<m;i+=1)
for(j=0;j<=i;j+=1)
	z=0
		for(k=0;k<n;k+=1)
			z+=xv[i][k]*xv[j][k]
		endfor
	xx[i][j]=z
	xx[j][i]=z
endfor
endfor
// X*Y'	
for(i=0;i<m;i+=1)
for(j=0;j<ny;j+=1)
	z=0
		for(k=0;k<n;k+=1)
			if(!numtype(yv[j][k]))
				z+=xv[i][k]*yv[j][k]
			endif
		endfor
	yx[i][j]=z
endfor
endfor
//
MatrixLLS /M=1 xx yx
duplicate /O M_B beta
MatrixTranspose beta
MatrixMultiply beta,xv
duplicate /O M_product yp
Variable /G rmsw=crms(yv,yp)
Variable /G r2w=CalcWilksFull(yv,yp) 
//return r2w // Wilks	
return rmsw
end



function CheckMatrixMultiply(a,b,s)
Wave a,b
String s
Variable n=DimSize(a,1)
Variable m=DimSize(b,0)
if(n!=m)
	print "matrix multiply mismatch at "+s
	Abort
endif
end

function CalcWilks(yv,yp,opt)  // Wilks 
Variable opt
Wave yv,yp
duplicate /O yp ep
ep=yv-yp
MatrixMultiply yv, yv /T
duplicate /O M_product yy
MatrixMultiply ep, ep /T
duplicate /O M_product ee
Variable r2=MatrixDet(ee)/MatrixDet(yy)// Wilks or Hotelling-Rozeboom
Variable n=DimSize(yy,0)
if(opt)
	Make /N=(n) /R /O r2vp
	r2vp=1-ee[p][p]/yy[p][p] // diagonals for n=1 give the usual R^2
endif
return r2
end

function CalcWilksFull(yv,yp)  // Wilks 
Wave yv,yp
Variable n=DimSize(yv,1)
Variable ny=DimSize(yv,0)
Make /N=(ny,ny) /O yy,ee
Make /N=(ny) /R /O r2vp
Variable i,j,k,z,zp		
for(i=0;i<ny;i+=1)
for(j=0;j<=i;j+=1)
	z=0
	zp=0
		for(k=0;k<n;k+=1)
			if( (!numtype(yv[i][k])) && (!numtype(yv[j][k])) )
				z+=yv[i][k]*yv[j][k]
				zp+=(yv[i][k]-yp[i][k])*(yv[j][k]-yp[j][k])
			endif
		endfor
	yy[i][j]=z
	yy[j][i]=z
	ee[i][j]=zp
	ee[j][i]=zp
endfor
endfor		
//
Variable r2=MatrixDet(ee)/MatrixDet(yy)
r2vp=1-ee[p][p]/yy[p][p] // diagonals for n=1 give the usual R^2
return r2
end

function crms(yv,yp)
Wave yv,yp
Variable ny=DimSize(yv,0)
Variable  n=DimSize(yv,1)
duplicate /O yv,ep
ep=yv-yp
Variable z=0,i,j,m=0
for(i=0;i<ny;i+=1)
for(j=0;j<n;j+=1)
		if(!numtype(ep[i][j]))
			z+=ep[i][j]^2
			m+=1
		endif
endfor
endfor
return sqrt(z/m)
end

//Window Table_test() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,42.5,510,218) yv,yp
	ModifyTable width(yv)=68,elements(yv)=(-3,-2),elements(yp)=(-3,-2)
EndMacro

//function SetYShowProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
Wave yv,yp,xv
Variable /G yshow
Variable m=DimSize(yv,1)
Variable n=DimSize(yv,0)
if(yshow<0)
	yshow=n-1
endif
if(yshow>=n)
	yshow=0
endif
Make /N=(m) /O /R yp1,yv1
Make /N=(DimSize(xv,0)) /O /R beta1,beta2
yv1=yv[yshow][p]
yp1=yp[yshow][p]
beta1=beta[yshow][p]
beta2=abs(beta1)
Variable /G r2show=r2vp[yshow]
End



Window show_y_plot() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,42.5,453,469.25) line45 vs line45
	AppendToGraph yv1 vs yp1
	ModifyGraph mode(yv1)=3
	ModifyGraph marker(yv1)=8
	ModifyGraph lSize(yv1)=3
	ModifyGraph rgb(line45)=(0,0,0)
	ModifyGraph msize(yv1)=2
	ModifyGraph opaque(yv1)=1
	ModifyGraph grid=2
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=25
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph lblLatPos(bottom)=-1
	Label left "exp."
	Label bottom "estimate"
	SetAxis left -4,4
	SetAxis bottom -4,4
	ShowInfo
	ShowTools
	SetVariable setvar0,pos={364,514},size={90,16},proc=SetYShowProc,title="show y"
	SetVariable setvar0,value= yshow
	SetVariable setvar1,pos={460,514},size={90,16},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
EndMacro




//================== logistic regression ==============================


// yv is binary
// xv is unshifted and it is (1,x1,x2,....) => beta0,beta1,beta2...

function LogisticReg_(xv,yv)
Wave xv,yv
Wave M_product,M_B
duplicate /O yv,pv,pp
duplicate /O xv,xw
Variable m=DimSize(xv,0)
Variable n=DimSize(xv,1)
Make /N=(1,m) /O beta
Variable /G lambda // =1e-3 defauly
if(!lambda)
	lambda=1e-3 
endif
Variable nt=100 // iterations
Variable k=0
beta=0
do
		MatrixMultiply beta,xv
		duplicate /O M_product yp
		pv=Logit(yp)
		xw=xv[p][q]*(pv[0][q]*(1-pv[0][q])+lambda*(p==q))
		pv=yv-pv
		MatrixMultiply xw, xv /T
		duplicate /O M_product xx
		MatrixMultiply xv, pv /T
		duplicate /O M_product yx
		MatrixLLS /M=1 xx yx
		MatrixTranspose M_B
		beta+=M_B
		if(Logistic_criterion(M_B,beta))
			//printf "converged after %d iteration\r\n",k
			return k
		endif
	k+=1
while(k<nt)
return -1
end

function LogisticRegB_(xvb,yvb)
Wave xvb,yvb
Wave beta
Wave M_product,M_B
MatrixMultiply beta,xvb
duplicate /O M_product ypb
Variable n=Logistic_n(yvb) 
Variable z0=Logistic_L0(yvb) 
Variable z=Logistic_qual(yvb,ypb)
Variable /G rmswb  = -2*z
z*=n/z0
Variable /G r2wb=z
end


function Logistic_r2() // McFadden pseudo-R2
Wave yv,yp
Variable n=Logistic_n(yv) 
Variable z0=Logistic_L0(yv) 
Variable z=Logistic_qual(yv,yp)
Variable /G rmsw = -2*z  // mean deviance
z*=n/z0
Variable /G r2w=z
Make /N=1 /O /R r2vp
r2vp[0]=1-z
return z
end

function Logistic_n(yv) 
Wave yv
Variable i,n=DimSize(yv,1)
Variable m=0
for(i=0;i<n;i+=1)
	if(!numtype(yv[0][i]))
		m+=1
	endif
endfor
return m
end

function Logistic_L0(yv) 
Wave yv
Variable i,n=DimSize(yv,1)
Variable n1=0,m=0
for(i=0;i<n;i+=1)
	if(!numtype(yv[0][i]))
		n1+=yv[0][i]
		m+=1
	endif
endfor
Variable n0=m-n1
return n0*ln(n0/m)+n1*ln(n1/m)
end

function Logistic_qual(yv,yp) 
Wave yv,yp
Variable i,n=DimSize(yv,1)
Variable z=0,m=0
for(i=0;i<n;i+=1)
	if(!numtype(yv[0][i]))
		z+=yv[0][i]*yp[0][i]-Logistic_Log(yp[0][i])	
		m+=1
	endif
endfor
return z/m
end

function Logit(z)
Variable z
if(z<-20)
	return 0
else
	return 1/(1+exp(-z))
endif
end


function Logistic_criterion(db,b)
Wave db,b
Wave yv
Variable criterion=1e-5
Variable m=DimSize(b,1)
Variable z1=0,z2=1,j
for(j=0;j<m;j+=1) // convergence of beta
	z1+=abs(db[0][j])
	z2+=abs(b[0][j])
endfor
if(z1/z2<criterion)
	return 1
else
	return 0
endif
end

function Logistic_dev(a,b)
Variable a,b
Variable d=a*Logistic_Log(-b)+(1-a)*Logistic_Log(b)
d=sqrt(abs(2*d))
if(a)
	return d
else
	return -d
endif
end

function Logistic_Log(u) 
Variable u
	if(u>10)
		return u
	endif
	if(u<-10)
		return exp(u)
	endif
return ln(1+exp(u))
end


// ======================== aide for Ridge regularization =============================


function hat_eigen_MLR()
Wave xv
MatrixMultiply xv, xv /T
MatrixEigenV /SYM /O M_product // X'X matrix, Hat matrix H = X (X'X)^-1 X'
duplicate /O W_eigenvalues hv
hv=abs(hv)
Variable /G lev=WaveMin(hv)
end

function hat_eigen_LGR()
Wave xv,xw,pv
xw=xv[p][q]*pv[0][q]*(1-pv[0][q])
MatrixMultiply xw, xv /T
MatrixEigenV /SYM /O M_product
duplicate /O W_eigenvalues hv
hv=abs(hv)
Variable /G lev=WaveMin(hv)
end


Window hat_matrix_ev() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(612,52,784,308) hv
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph msize=2
	ModifyGraph mirror=2
	ModifyGraph minor(left)=1
	ModifyGraph axOffset(left)=-0.6
	Label left "eigenvalue (hat matrix)"
	Label bottom " "
	SetAxis/A/N=1/E=1 left
	SetAxis/A/N=1/E=1 bottom
	Cursor/P A hv 5
	ShowInfo
EndMacro



