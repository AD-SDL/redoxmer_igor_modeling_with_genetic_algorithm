// coefficients validation : partitioning w/o parameter space optimization


Proc CoeffValidation(n,part)
Variable n=100,part=80
Prompt n,"number of trials"
Prompt part,"train partition, say 70%"
Silent 1
PauseUpdate
Variable /G do_regression=1
Variable /G rmsw,rmswb,r2w,r2wb
Variable i=0,a,b
Make /N=(n) /O /R r2a,r2b,r1a,r1b
Variable m=DimSize(yval_set,1)
r2a=nan;r2b=nan;r1a=nan;r1b=nan
part/=100
//
duplicate /O status_set status_set_old
choose_y(yname_set)
binary_weight_fixed(1)
EvVect_(VB,0)
Variable /G r2full=1-r2w
Variable /G rmsfull=rmsw
printf ">> full set r2=%.3f rms=%.4f\r\n\r\n",r2full,rmsfull
DoWindow /K cross_val_rms
DoWindow /K cross_val_R2
//
do
	binary_weight_fixed(part)
	EvVect_(VB,0)       // training
	b=1-EvVect_(VB,1)   // validation
	r2a[i]=1-r2w
	r2b[i]=1-r2wb
	r1a[i]=rmsw
	r1b[i]=rmswb
	printf "trial %d: r2(train)=%.3f r2(val)=%.3f; rms(train)=%.4f rms(val)=%.4f\r\n",1+i,1-r2w,1-r2wb,rmsw,rmswb
	if(ButtonIsClicked())
		break
	endif
i+=1
while(i<n)
//
duplicate /O status_set_old status_set
binary_weight_fixed(1)
EvVect_(VB,0)
//
printf "\r\n>>> full set R2=%.3f :: %d:%d partitioning\r\n",r2full,100*part,100*(1-part)
printf ">>> medians validation: R2=%.3f, rms=%.4f\r\n",r2b[n/2],r1b[n/2]
Variable ava,stda,avb,stdb
WaveStats /Q r2a
ava=V_avg; stda=V_sdev
WaveStats /Q r2b
avb=V_avg; stdb=V_sdev
printf ">>> R2(train)=%.3f±%.3f; R2(validation)=%.3f±%.3f\r\n",ava,stda,avb,stdb
WaveStats /Q r1a
ava=V_avg; stda=V_sdev
WaveStats /Q r1b
avb=V_avg; stdb=V_sdev
printf ">>> rms(train)=%.3f±%.3f; rms(validation)=%.3f±%.3f\r\n",ava,stda,avb,stdb
//
Sort r2b,r2a,r2b
Sort r1b,r1a,r1b
r1a/=rmsfull
r1b/=rmsfull
Make /N=(n/10) /O /R d2a,d2b,d1a,d1b
Histogram/B=1 r2a,d2a
Histogram/B=1 r2b,d2b
Histogram/B=1 r1a,d1a
Histogram/B=1 r1b,d1b
coef_val_rms()
coef_val_R2()
end

Proc CVDissmissButton(ctrlName) : ButtonControl
	String ctrlName
DoWindow /K coef_val_rms
DoWindow /K coef_val_R2
End


//================== cross validation by partition ======================================//


Proc CrossValidation(n,part,NTmax,opt)
Variable n=100,part=80,NTmax=1500,opt=0
Prompt n,"number of trials"
Prompt part,"train partition, say 70%"
Prompt NTmax,"maximum iteration"
Prompt opt,"metric: 0 -mine, 1 - Hausdorff, 2 - Binet-Cauchy"
Silent 1
PauseUpdate
//
Variable /G VB_metric=opt
Variable /G beta_cons_offset=0
Variable /G reg_type
Variable /G Ng=numpnts(VB)
Variable /G NT=NTmax
Variable /G rmsw,r2w,rmswb,r2wb
Variable i
Variable nhist=20
//
duplicate /O VB VB_old
duplicate /O status_set status_set_old
//
Make /N=0 /O /R VB_score,VB_score_corr,VB_wt
Make /N=0 /O /T VB_valid,VB_cons_comment
Make /N=0 /O /R VB_cons,beta_cons,rms_t_cons,rms_v_cons,r2_t_cons,r2_cons
Make /N=0 /O /R VB_cons_all,beta_cons_all
Make /N=(Ng+1) /O /R score_valid
Make /N=(nhist) /O /R score_corr_valid
SetScale/I x 0,1,"", score_corr_valid,score_valid
score_valid=0
score_corr_valid=0
DoWindow /K mismatch_count_validation
mismatch_count_validation()
//
part/=100
choose_y(yname_set)
//
do	
	
	printf "\r\n\r\n =========== partition %d =========== \r\n\r\n",i+1
	
	if(i)
		Prepare_Validation(part)
		Genetic(NT)	
	endif
	
	if(consensus_input())
		printf "\r\nRMS (training)=%.4f, (validation)=%.4f :: R2 (training) = %.4f\r\n",rmsw,rmswb,1-r2w
	else
		printf "\r\nalready logged"
	endif
	
	if(ButtonIsClicked())
		break
	endif
i+=1
while(i<n)
//
FinishCrossValidation(i+1)
status_set=status_set_old
Prepare_Validation(1)
subset_train_valid()
duplicate /O status_set color_set
Variable /G solution_found=1
ShowBestSolution()
r2_cons=1-r2w
DoWindow /K Table_coeff
end


Proc FinishCrossValidation(i)
Variable i
Variable /G reg_type
Variable m=DimSize(beta,0)
Variable n=DimSize(beta,1)-reg_type
Variable /G beta_cons_offset
//
VB_wt/=i
score_valid/=i
score_corr_valid/=i
beta_cons/=i
beta_cons_offset/=i
//
Redimension /N=(m,numpnts(VB_cons)) beta_cons
sort_beta_cons()
//
if(reg_type) 
   m=numpnts(VB_wt)
	Redimension /N=(n+1,m) beta_cons_all
	Redimension /N=(n,m)   VB_cons_all
endif
end


Proc Prepare_Validation(part)
Variable part
	Variable /G r2show=nan
	Variable /G solution_found=0
	Variable /G NT
	binary_weight_fixed(part)
	choose_fixed_subset()
	duplicate /O status_set color_set
	VB=VB_old
	Prepare_Genetics(NT)		
end

function mismatch_valid()
Wave VB,VB_old
Variable i,j,n=numpnts(VB)
Variable score=n
for(i=0;i<n;i+=1)
	for(j=0;j<n;j+=1)
		if(VB[i]==VB_old[j])
			score-=1
			break
		endif
	endfor
endfor
return score/n // 0 = perfect match, 1 - worst match
end


function dval_d_max(i)
Variable i
Wave CM
Variable ro=0,j,n=DimSize(CM,0)
for(j=0;j<n;j+=1)
	ro=max(ro,abs(CM[i][j]))
endfor
return 1-ro
end


function mismatch_corr_valid(opt)
Variable opt
Wave VB,VB_old
Variable i,n=numpnts(VB)
Variable score=0
Make /N=0 /O /R M_B
Make /N=(n,n) /O /R CA,CB,CM
CM=Correlation_2D(dval_set,VB[p],VB_old[q],1)
CA=Correlation_2D(dval_set,VB[p],VB[q],1)
CB=Correlation_2D(dval_set,VB_old[p],VB_old[q],1)
//
switch(opt)
	case 0: // my score 
		for(i=0;i<n;i+=1)
			score+=dval_d_max(i)/n
		endfor
	break

	case 1: // Hausdorff 
	for(i=0;i<n;i+=1)
		score=max(score,dval_d_max(i))
	endfor
	break

	case 2: //Binet-Cauchy 	 
		score=1-abs(MatrixDet(CM))/sqrt(MatrixDet(CA)*MatrixDet(CB))
	break
	
	case 3: //projection	
	default: 
		score=1-norm(CM)^2/(norm(CA)*norm(CB))
	break
	
endswitch

return score

end


function /S VB_string()
Wave VB
Variable i,n=numpnts(VB)
String s=""
for(i=0;i<n;i+=1)
	s+=num2str(VB[i])+";"
endfor
return SortList(s,";",2)
end

function consensus_param_input()
Variable /G VB_metric
Wave score_valid,score_corr_valid
Wave /T VB_valid
Wave VB_count,VB_score,VB_score_corr,VB_wt
Wave rms_t_cons,rms_v_cons,r2_t_cons,r2_cons
//
Variable /G rmsw,rmswb,r2w
Variable j,i=numpnts(VB_valid)
String s = VB_string()
Variable k=WhereInSWave(VB_valid,s)

//
if(k<0)
	Make /N=(i+1) /O /R VB_count,VB_score,VB_score_corr,VB_wt
	Make /N=(i+1) /O /R rms_t_cons,rms_v_cons,r2_t_cons,r2_cons
	Make /N=(i+1) /O /T VB_valid
	
	VB_valid[i]=s
	VB_wt[i]=1
	
	j=mismatch_valid()
	score_valid[x2pnt(score_valid,j)]+=1	
	VB_score[i]=j	
	
	j=mismatch_corr_valid(VB_metric)
	score_corr_valid[x2pnt(score_corr_valid,j)]+=1	
	VB_score_corr[i]=j	
	
	rms_t_cons[i]=rmsw
	rms_v_cons[i]=rmswb	
	r2_t_cons[i]=1-r2w
	
	return 1
else

	VB_wt[k]+=1
	
	j=x2pnt(score_valid,VB_score[k])
	score_valid[j]+=1	

	j=x2pnt(score_corr_valid,VB_score_corr[k])
	score_corr_valid[j]+=1
	
	return 0	
endif  	
end

function consensus_input()
Variable /G reg_type
Variable /G beta_cons_offset
Wave /T dname_set,VB_cons_comment
Wave VB,VB_cons,beta,beta_cons
Wave VB_cons_all,beta_cons_all
//
Variable i,j,k,j0,n0,flag
Variable m=Dimsize(beta,0)
Variable n=Dimsize(beta,1)-reg_type
Variable nc=numpnts(VB_cons)
//
if(reg_type)
	beta_cons_offset+=beta[0][0]
endif
for(i=0;i<n;i+=1)
	k=WhereInWave(VB_cons,VB[i])
	if(k<0)
		j0=nc*m		
		Make /N=(nc+1) /O /R VB_cons
		Make /N=(nc+1) /O /T VB_cons_comment
		Make /N=(j0+m) /O /R beta_cons
		VB_cons[nc]=VB[i]
		VB_cons_comment[nc]=dname_set[VB[i]]
		for(j=0;j<m;j+=1)
			beta_cons[j0+j]=beta[j][i+reg_type]
		endfor
		nc+=1
	else
		j0=k*m
		for(j=0;j<m;j+=1)
			beta_cons[j0+j]+=beta[j][i+reg_type]
		endfor
	endif
endfor
//
flag=consensus_param_input()
if(reg_type && flag)	
	n0=numpnts(beta_cons_all)
	Make /N=(n0+n+1) /O /R beta_cons_all
	beta_cons_all[n0,n0+n]=beta[0][p-n0]
	n0=numpnts(VB_cons_all)
	Make /N=(n0+n) /O /R VB_cons_all
	VB_cons_all[n0,n0+n-1]=VB[p-n0]
endif	
return flag
end


function sort_beta_cons()
Wave beta_cons,VB_cons
Wave /T VB_cons_comment
variable /G Ng
Variable i,j,z
Variable m=DimSize(beta_cons,0)
Variable n=DimSize(beta_cons,1)
Make /O /N=(n) /R beta_abs_cons
Make /O /N=(n) /I beta_cons_indx,indx
beta_cons_indx=0
beta_cons_indx[0,Ng-1]=1
for(i=0;i<n;i+=1)
	z=0
	for(j=0;j<m;j+=1)
		z+=beta_cons[j][i]^2
	endfor
	beta_abs_cons[i]=sqrt(z/m)
endfor
MakeIndex /R beta_abs_cons,indx
IndexSort indx,beta_abs_cons,VB_cons,VB_cons_comment,beta_cons_indx
duplicate /O beta_cons beta_cons_
beta_cons_=beta_cons[p][indx[q]]
duplicate /O beta_cons_ beta_cons
KillWaves beta_cons_,indx
end


Proc RestoreOriginalConsensus()
restore_from_cons()
BestVectorFromList()
end

function restore_from_cons()
Wave beta_cons_indx
Wave /T VB_cons_comment
Variable i,m=0,n=numpnts(beta_cons_indx)
for(i=0;i<n;i+=1)
	if(beta_cons_indx[i])
		Make /N=(m+1) /O /T V_comment
		Make /N=(m+1) /O /R VB
		V_comment[m]=VB_cons_comment[i]
		m+=1
	endif
endfor
end


Proc ImportConsensus(n)
Variable n=-1
Prompt n,"import n highest consensus variable; -1 for all"
Variable /G reg_type
Variable /G solution_found=1
if(!reg_type)
	consensus_solution_MLR(n)
	SetYShowProc_(1)	
	DoWindow /K Table_coeff
	Table_coeff()
else
	consensus_solution_LGR()
endif					
color_comp_set()
dexplain_show=des_explain(V_comment_show)
end


function consensus_solution_MLR(n)
Variable n
Wave VB_cons,beta_cons
Wave /T dname_set
Variable m=DimSize(beta_cons,0)
if(n>0)
	Make /O /N=(n) /R VB
	Make /N=(m,n) /O /R beta
	VB=VB_cons[p]
	beta=beta_cons[p][q]
else
	duplicate /O VB_cons VB
	duplicate /O beta_cons beta
endif
Variable /G Ng=numpnts(VB)
Make /N=(Ng) /O /T V_comment
V_comment=dname_set[VB]
duplicate /O V_comment,V_comment_show
Wave /I subset_train,subset_valid
Make /N=0 /O /R xv,yv,xvb,yvb,yp,ypb
//
			if(AutoFill(VB,xv,subset_train,0))
				AutoShiftY(yv,ysubset,subset_train,0)
				MatrixMultiply beta,xv
				duplicate /O M_product yp
			endif
				

			if(AutoFill(VB,xvb,subset_valid,1))
				AutoShiftY(yvb,ysubset,subset_valid,1)
				MatrixMultiply beta,xvb
				duplicate /O M_product ypb
			endif
			
end


function consensus_solution_LGR()
Wave VB_cons_all,beta_cons_all,VB_wt
Wave yval_set
Wave /T yname_set,dname_set
Wave /I subset_train,subset_valid,ysubset
Make /N=0 /O /R xv,yv,xvb,yvb,yp,ypb,pv
//
Variable /G yshow=0
String /G yname_show=yname_set[0]
Variable /G Ng=DimSize(VB_cons_all,0)
Variable m=numpnts(VB_wt)
Make /N=(Ng) /O /R VB
Make /N=(Ng) /T /O V_comment
Make /N=(1,Ng+1) /O /R beta
Variable i
//
Variable l=DimSize(yval_set,1)
Make /N=(l) /O /R ep1,yp1,yv1,ys1,tv1,pv1
ys1=in_train_set(p)	
yv1=yval_set[0][p]
pv1=0
ep1=nan
//		
for(i=0;i<m;i+=1)

	VB=VB_cons_all[p][i]
	V_comment=dname_set[VB]
	beta[0][]=beta_cons_all[q][i]
	
			if(AutoFill(VB,xv,subset_train,0))
				AutoShiftY(yv,ysubset,subset_train,0)
				MatrixMultiply beta,xv
				duplicate /O M_product yp
			endif
				
			if(AutoFill(VB,xvb,subset_valid,1))
				AutoShiftY(yvb,ysubset,subset_valid,1)
				MatrixMultiply beta,xvb
				duplicate /O M_product ypb
			endif
					
		fill_yp1(0)	
	   pv1+=VB_wt[i]*Logit(yp1)
			
endfor
//
tv1=label_text(yv1,pv1)
VB=VB_cons_all[p][0]
V_comment=dname_set[VB]
beta[0][]=beta_cons_all[q][0]
Make /O /N=0 /R beta1,beta2
Variable /G r2show=nan
End

Macro test()
	consensus_input()
	FinishCrossValidation(0)
end

//=========================================================================//


Window coef_val_rms() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(73,49,478,491) d1b,d1a
	ModifyGraph mode=5
	ModifyGraph rgb(d1b)=(0,0,52224)
	ModifyGraph hbFill(d1a)=2
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axOffset(bottom)=0.458333
	ModifyGraph axThick=2
	Label left "frequency"
	Label bottom "rms/rms(all)"
	SetAxis/A/N=1/E=1 left
	SetAxis/A/N=1 bottom
	ShowInfo
	Legend/C/N=text0/J/F=0/B=1/A=MC/X=31.53/Y=42.32 "\\Z08sets:\r\\s(d1a) training\r\\s(d1b) validation"
EndMacro

Window coef_val_R2() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(489,46,903,491) d2a
	ModifyGraph mode=5
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axOffset(bottom)=0.458333
	ModifyGraph axThick=2
	Label left "frequency"
	Label bottom "training R\\S2"
	SetAxis/A/N=1/E=1 left
	SetAxis/N=1 bottom*,1
	Cursor/P A d2a 66
	ShowInfo
	Button button0,pos={343.00,417.00},size={50.00,20.00},proc=CVDissmissButton,title="dismiss"
EndMacro




Window rms_consensus() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(660,502,999,806) rms_v_cons vs rms_t_cons
	ModifyGraph mode=3
	ModifyGraph marker=8
	ModifyGraph lSize=2
	ModifyGraph msize=3
	ModifyGraph opaque=1
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axThick=2
	Label left "rms (validation)"
	Label bottom "rms (training)"
	SetAxis/A/N=1 left
	SetAxis/A/N=1 bottom
	ShowInfo
EndMacro

Window r2_consensus() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(660,502,999,806) r2_cons,r2_t_cons
	ModifyGraph mode(r2_t_cons)=3
	ModifyGraph marker(r2_cons)=8,marker(r2_t_cons)=19
	ModifyGraph lSize=2
	ModifyGraph lStyle(r2_cons)=3
	ModifyGraph rgb(r2_cons)=(0,0,0)
	ModifyGraph msize=3
	ModifyGraph opaque=1
	ModifyGraph mirror=2
	ModifyGraph minor(left)=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axThick=2
	Label left "R\\S2"
	SetAxis/A/N=1 left
	SetAxis/A/N=1 bottom
	ShowInfo
EndMacro


Window coef_consensus() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(393,54.5,897.75,552.5) beta_abs_cons
	ModifyGraph mode=1
	ModifyGraph lSize=2
	ModifyGraph zColor(beta_abs_cons)={beta_cons_indx,*,*,Rainbow,1}
	ModifyGraph log(left)=1
	ModifyGraph mirror=2
	ModifyGraph minor(bottom)=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axThick=2
	Label left "absolute MLR coefficients"
	Label bottom "descriptors"
	SetAxis left 0.001,*
	ShowInfo
	TextBox/C/N=text0/F=0/A=MC/X=28.47/Y=37.41 "\\Z09consensus coefficients\rcross validation\rred - whole data\r(model descriptors)\rpurple - leave 10% off"
EndMacro


Window Table_consensus_score() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(481,340,1165,713) VB_wt,VB_score,VB_score_corr,VB_valid,beta_cons_all
	ModifyTable format(Point)=1,width(Point)=0,alignment(VB_wt)=0,format(VB_wt)=3,width(VB_wt)=48
	ModifyTable alignment(VB_score)=0,format(VB_score)=3,width(VB_score)=50,rgb(VB_score)=(39321,1,1)
	ModifyTable alignment(VB_score_corr)=0,format(VB_score_corr)=3,width(VB_score_corr)=44
	ModifyTable alignment(VB_valid)=0,width(VB_valid)=130,rgb(VB_valid)=(1,4,52428)
	ModifyTable elements(beta_cons_all)=(-3,-2)
EndMacro


Window Table_consensus_beta() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(811,258,1222,879) VB_cons,VB_cons_comment,beta_cons
	ModifyTable format(Point)=1,width(Point)=28,alignment(VB_cons)=0,width(VB_cons)=52
	ModifyTable rgb(VB_cons)=(0,0,65535),alignment(VB_cons_comment)=0,width(VB_cons_comment)=94
	ModifyTable rgb(VB_cons_comment)=(65535,0,0),format(beta_cons)=3,width(beta_cons)=54
	ModifyTable elements(beta_cons)=(-3,-2)
EndMacro


Window mismatch_count_validation() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(612,52,784,308) score_corr_valid,score_valid
	ModifyGraph mode(score_corr_valid)=5,mode(score_valid)=1
	ModifyGraph lSize(score_valid)=5
	ModifyGraph rgb(score_corr_valid)=(1,4,52428)
	ModifyGraph mirror=2
	ModifyGraph axOffset(left)=-2.4
	Label bottom "mismatch score"
	SetAxis/A/N=1/E=1 left
EndMacro
