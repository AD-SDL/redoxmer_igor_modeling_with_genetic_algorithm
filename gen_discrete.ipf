// (c) Eli Shkrob Igorized 3/25/03

// initial guess and best vector is VB
// dV 	is variablily per generation; if 0 - no variation; if <0 strong mutations
//		can reverse the sign

//	Nv 		number of "vectors"
//	Ng		number of "genes"
//	V		matrix of Nv vectors
//	VN		new matrix
//	Nvs		=> Nvs/Nv = survival rate
//	C		fitness
//	Nst		number of mutations per vector
//	VC		current vector
//	VB		the fittest vector
// 	NT 		is the number of iterations
// rmsB	is the best rms

//Macro DoTestGenetics()
Silent 1
	Variable /G DBG=0
	Variable /G Nv=50
	Variable /G Nvs=10
	Variable /G Nst=2
	Variable /G NT=500
	Variable /G rmsB=0
	Variable /G IT=0
	Genetic(NT)
	END

//		genetic algorithm for optimization
//           discrete number version 3-6-2018

function Genetic(n)
Variable n
Variable /G NT=n
Variable /G DBG=0		// debugging
Wave VC,VB,VN,V
Variable /G Nv,Ng,Nvs,Nst,IT=1
Variable /G rmsB=1e20
Variable i
Wave subset
Variable /G Ng
Make /O /N=(Ng) /R VC,VB
Make /O /N=(Nv,Ng) /R VN,V
Make /O /N=(Nv) /R C,PW,IC
Make /T /O /N=(Ng) V_comment
// random vectors
	for(i=0;i<Nv;i+=1)
		choose_random_subset()
		V[i][]=subset[q]
	endfor
	GeneticContinue(NT)
end
	
function GeneticContinue(NT)
Variable NT
Variable /G DBG,rmsB,IT=1
//
Variable i
Wave V,VB,VC
V[0][]=VB[q]  // keep the best vector
for(i=1;i<=NT;i+=1)
	IF(DBG)
		printf "\r\r>>> GENERATION  %d :: best rms=%e",IT,rmsB
	endif	
	if(OneIteration())
		break
	endif
	IT+=1	
	if(ButtonIsClicked())
		break
	endif
endfor
	VC=VB
	EvVect(VB)
	renew_best()
END

function OneIteration()
Variable i,j,z
Variable /G Nst,rmsB,Ng,Nv,Nvs
Variable /G DBG
Variable M,N,MC,NC,IX
Wave     VC,VB,VN,V,C,IC,PW
Variable mutations=0
//   mutation

	IF(DBG)
		printf "\raverage of %d mutations per vector of %d genes\r\n",Nst,Ng
		printf "\rnew vectors evaluated:"
	endif
	
	i=0
	for(i=0;i<Nv;i+=1)
		VC=V[i][p]	
		if(i)
			mutations=point_mutation(VC,Nst)
		endif
		IF(DBG) 
		printf "\r\n%d mutations out of average of %d",mutations,Nst
		printf "\r\nvector=%d :: ",i
			for(j=0;j<Ng;j+=1)
				printf "%d," VC[j]
			endfor
		endif
		z=EvVect(VC)
		V[i][]=VC[q]
		if(DBG)
			printf ":: rms=%g",z			
		endif	
		if(z<rmsB) 
			rmsB=z
			VB=VC
			renew_best()
			doUpdate
		endif
	
	C[i]=z
	i+=1
	endfor

	IC=p
	Sort C,C,IC
	
	IF(DBG) 
		print "\rsorted chi2 array :: number, element, index"
		i=0
		do
			printf "\r%d %g %d", I,C[i],IC[i]
			i+=1
		while(i<Nv)
		printf "\rrecombination cross-table:"
	endif

		PW=1/C
		if(numtype(PW[0]))
			return 1
		endif
		if(PW[0]>1e10)
			return 1
		endif

for(i=0;i<Nv;i+=1)
		M=IPL(Nvs,PW)
		do
			N=IPL(Nvs,PW)
		while(N==M) 
		MC=IC[M]
		NC=IC[N]
		IX=irand(Ng)
	
	IF(DBG) 
		printf "\rvector %d <= %d + %d; crossed @ %d",I,MC,NC,IX
	endif
	
//	recombination 
	      VN[i][0,IX]=V[MC][q]
		VN[i][IX+1,Ng-1]=V[NC][q]
		
endfor

		V=VN
		return 0
end
	

//	pick i for a discrete distribution P
	
	FUNCTION IPL(N,P)
	Wave P
	Variable N
	Variable z=rand()*sum(P,0,N-1)
	Variable x=0,y=P[0],IPL=0

	do
	IF(((z>x) %& (z<y)) %| (IPL==N-1))
		return IPL
	endif
	IPL+=1
	x=y
	y+=P[IPL]
	while(1)
END
	
Function rand()
	return .5*(1+enoise(1))
end

Function irand(N)
Variable N
	return floor(N*rand())
end

Window Table_vector() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(797,48,1019,346) VB,V_comment
	ModifyTable format(Point)=1,width(Point)=26,alignment(VB)=0,format(VB)=1,digits(VB)=2
	ModifyTable width(VB)=48,rgb(VB)=(1,4,52428),alignment(V_comment)=0,width(V_comment)=132
	ModifyTable rgb(V_comment)=(65535,0,0)
EndMacro





Window GeneticsPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1002,748,1227,954) as "Genetic Algorithm Optimization"
	SetVariable Nv_control,pos={17.00,8.00},size={150.00,18.00},title="# of vectors"
	SetVariable Nv_control,font="Arial",fSize=12,format="%d"
	SetVariable Nv_control,limits={1,1000,1},value= Nv
	SetVariable Nst_control,pos={17.00,54.00},size={180.00,18.00},title="# strong mutations"
	SetVariable Nst_control,font="Arial",fSize=12,format="%d"
	SetVariable Nst_control,limits={1,1000,1},value= Nst
	SetVariable NVs_control,pos={17.00,77.00},size={180.00,18.00},title="# survived/g-n"
	SetVariable NVs_control,font="Arial",fSize=12,format="%d"
	SetVariable NVs_control,limits={1,1000,1},value= Nvs
	SetVariable NT_control,pos={17.00,100.00},size={180.00,18.00},title="# of iterations"
	SetVariable NT_control,font="Arial",fSize=12,format="%g"
	SetVariable NT_control,limits={5,100000,5},value= NT
	ValDisplay valdisp0,pos={17.00,32.00},size={150.00,17.00},title="# of parameters"
	ValDisplay valdisp0,font="Arial",fSize=12,format="%d",frame=0
	ValDisplay valdisp0,limits={0,0,0},barmisc={0,1000},value= #"Ng"
	Button go_button,pos={20.00,126.00},size={50.00,20.00},proc=DoGenetics,title="Go!"
	Button cont_button,pos={19.00,152.00},size={55.00,20.00},proc=ContinueGenetics,title="Continue"
	ValDisplay it_control,pos={86.00,128.00},size={110.00,17.00},title="iteration"
	ValDisplay it_control,font="Arial",fSize=12,frame=0
	ValDisplay it_control,limits={0,0,0},barmisc={0,4000},value= #"IT"
	ValDisplay rms_control,pos={93.00,153.00},size={110.00,17.00},title="rms"
	ValDisplay rms_control,font="Arial",fSize=12,format="%.2e",frame=0
	ValDisplay rms_control,limits={0,0,0},barmisc={0,1000},value= #"rmsB"
	CheckBox check0,pos={19.00,182.00},size={73.00,17.00},title="auto quad"
	CheckBox check0,variable= auto_quad
	Button button0,pos={96.00,180.00},size={50.00,20.00},proc=EvOnce,title="Evaluate"
	CheckBox check1,pos={152.00,182.00},size={48.00,17.00},title="do LR",fSize=11
	CheckBox check1,variable= do_regression
EndMacro

Proc EvOnce(ctrlName) : ButtonControl
	String ctrlName
EvVect_(VB,0)
renew_best()
End


function  ButtonIsClicked()
	GetMouse
	if (V_flag & 2)
		return(1)
	else
		return(0)
	endif
end

Function ContinueGenetics(ctrlName) : ButtonControl
	String ctrlName
	Variable /G NT
	Variable /G do_regression=1
	GeneticContinue(NT)
End

Function DoGenetics(ctrlName) : ButtonControl
	String ctrlName
	Variable /G NT
	Variable /G do_regression=1
	Genetic(NT)
End
