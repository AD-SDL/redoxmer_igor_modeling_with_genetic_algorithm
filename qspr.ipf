// set data as your folder and put all your descriptors into a subfolder ./des
// input files for matrix makers consist of compounds + smiles (CSV format all)
// 3-5-18

// * include in the set
// # exclude from the set

Menu "QSPR"
	Submenu "CVS import"
		"LoadCSV_des"
		"LoadCSV_table"
		"FindCorrespondence"
		"yCalculator"
		"yDelete"
		"yBinary"
	end
	Submenu "Descriptors"
		"CreateSetByList_des"
		"CreateSetAll_des"
		"CreateSet_HDF5_des"
		"ParseSet_des"
		"Test_salts"
		"Export_salts"
		"MkUniqueCompList"
		"ExportCompList"
		"ReplaceNansInDes"
	end
	Submenu "Hash and star"
		"ClassifyDescriptors"
		"HashSelection"
		"UnhashSelection"
		"HashAll2n3D"
		"HashAll3D"
		"HashSmallStd"
		"UnhashAllDescriptors"
		"UnstarAllDescriptors"
		"UnstarAllYs"
		"StarSelection"
		"StarAllYs"
		"SelectAllCompounds"
		"SelectAllNotNan"
	end
	Submenu "Sets"		
		"CopyOriginalSet"
		"RestoreOriginalSet"
		"RestoreLinSet"
		"RestoreQuadSet"
		"RestoreFinalSet"
		"AcceptSet"
		"AcceptLinear"
		"AcceptQuad"
		"AcceptFinal"
		"QuadExpansion"
		"SaveSet"
		"ExportSet"
	end
	Submenu "Regression"
		"MLR_Genetics"	
		"LGR_Genetics"
		"ExportL_LGR_result"
		"CoeffValidation"
		"TestDim"
		"TestDimLQ"
		"ErrorNeighborSort"
	end
	Submenu "Validation"
		"CrossValidation"
		"ImportConsensus"
		"RestoreOriginalConsensus"
	end
	Submenu "Evaluate"
				"EvalLoadedDesDir"
				"EvalLoadedSet"
				"EvalParamPlot"
	end
	"ShowBestSolution"
	"BestVectorFromList"
	"DoAllYPlots"
	"PlotAllMLR_Y"
	"MovieAllMLR_Y"
	"Copy_des_wave"
	"Hist_des_wave"
	"Copy_y_wave"
	"Hist_y_wave"
end

// all variables need to have names with ".y" and be at the end of the table
// alternatively load a separate (larger) list; correspondences will be tested

// assumes unshifted !!!


Window MainPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(612,343,783,460)
	ModifyPanel cbRGB=(65535,54607,32768)
	Button button0,pos={6.00,6.00},size={50.00,20.00},proc=MMLR_Button,title="MMLR"
	SetVariable setvar0,pos={70.00,6.00},size={50.00,19.00},title="n",format="%d"
	SetVariable setvar0,limits={0,inf,1},value= Ng
	Button button1,pos={8.00,30.00},size={50.00,20.00},proc=LGR_Button,title="LGR"
	SetVariable setvar1,pos={60.00,31.00},size={110.00,19.00},title="lambda"
	SetVariable setvar1,format="%.1e",limits={0,inf,0.001},value= lambda
	Button button2,pos={33.00,56.00},size={100.00,20.00},proc=BestSolution_Proc,title="Best Solution"
	SetVariable setvar2,pos={20.00,82.00},size={130.00,19.00},disable=2,title="smallest ev"
	SetVariable setvar2,format="%.1e",limits={0,inf,1},value= lev,noedit= 1
EndMacro



Proc MMLR_Button(ctrlName) : ButtonControl
	String ctrlName
Variable /G Ng,lambda
MLR_Genetics(Ng,lambda)
End

Proc LGR_Button(ctrlName) : ButtonControl
	String ctrlName
Variable /G Ng,lambda
LGR_Genetics(Ng,lambda)
End

Proc BestSolution_Proc(ctrlName) : ButtonControl
	String ctrlName
Variable state=WinType("Table_coeff")
if(!state)
	ShowBestSolution()
else
	DoWindow /K Table_coeff
endif
End

Proc SelectAllCompounds()
status_set=1
duplicate /O status_set color_set,ys1
end

Proc SelectAllNotNan()
status_set=(numtype(yval_set[0][p])) ? 0 : 1
end


Proc LoadCSV_des()
Silent 1
PauseUpdate
LoadWave/J/M/A=inp/P=home/K=0 
String /G y_transform=""
String /G s_path,s_filename
String /G lib=removefileextension(s_filename)
Variable n=DimSize(inp0,0)-1
Variable nd=DimSize(inp0,1)
Make /N=(nd) /T /O kdes
kdes=inp0[0][p]
Variable /G ny=CountInSWave(kdes,".y")
nd-=2+ny
Variable ns=WhereInSWaveLoose(kdes,"status")
if(ns>=0)
	nd-=1
endif
printf "%d compounds x %d descriptors/parameters\r\n",n,nd
Make /N=(n) /O /T comp_set,smi_set
Make /N=(nd) /O /T dname_set
Make /N=(nd,n) /O /R dval_set
comp_set=inp0[p+1][0]
smi_set=inp0[p+1][1]
dname_set=inp0[0][p+2]
dval_set=str2num(inp0[q+1][p+2])
if(ny)
	printf "loading %d variables\r\n",ny
	Make /N=(ny,n) /O /R yval_set
	Make /N=(ny,ny) /O /R yval_ro
	Make /N=(ny) /O /T yname_set
	Make /N=(ny) /O /R r2y
	r2y=nan
	Make /N=(n) /O /B status_set
	if(ns<0)
		status_set=1
	else	
		status_set=str2num(inp0[p+1][ns])
	endif
	yname_set=inp0[0][p+2+nd]
	yval_set=str2num(inp0[q+1][p+2+nd])
	yval_ro=Correlation_2D(yval_set,p,q,1)
else
	Make /N=0 /O /R yval_set,r2y
	Make /N=0 /O /T yname_set
	Make /N=0 /O /B status_set
	DoAlert 0,"load the table of variables and parameters separately"
endif
KillWaves inp0,kdes
Variable /G table_inp=0
mkucomp()
Variable /G auto_quad=0
end

function mkucomp()
Wave /T comp_set
String /G ucomp_set=""
Variable i,m=0,n=numpnts(comp_set)
String s
for(i=0;i<n;i+=1)
	s=comp_set[i]
	if(WhichListItem(s,ucomp_set)<0)
		ucomp_set+=s+";"
		m+=1
	endif
endfor
printf ">>> %d entries :: %d compounds\r\n",n,m
end

function mkparam()
Wave /T dname_set
String /G par_set=""
Variable i,m=0,n=numpnts(dname_set)
String s
for(i=0;i<n;i+=1)
	s=dname_set[i]
	if(strsearch(s,"*",0)>=0)
		par_set+=s+";"
		m+=1
	endif
endfor
printf ">>> %d entries :: %d starred parameters\r\n",n,m
end


// Always load descriptors before loading  the table of parameters and variables (y)
// variables all have ".y" in their names

// assumes unshifted !!!

Proc LoadCSV_table(opt)
Variable opt=1
Prompt opt,"1 - keep missing compounds"
Silent 1
PauseUpdate
DoAlert 0, "load descriptors first"
LoadWave/J/M/A=inp/P=home/K=0 
String /G y_transform=""
Variable n=DimSize(inp0,0)-1
Variable /G auto_quad=0
if(!n)
	abort "no compounds"
endif
Variable m=DimSize(inp0,1)-1
Make /N=(m) /O /T kdes
kdes=inp0[0][p+1]
Variable  /G ny=CountInSWave(kdes,".y")
if(!ny)
	abort "no y-variables"
endif
Variable ns=WhereInSWaveLoose(kdes,"status")
Variable /G  nx=m-ny
if(ns>=0)
	nx-=1
endif
printf "%d compounds x %d  variables + %d parameters\r\n",n,ny,nx
//
Make /N=(n) /O /B status_set_inp
Make /N=(n) /O /T ycomp_set_inp
Make /N=(ny) /O /T yname_set
Make /N=(ny) /O /R r2y
r2y=nan
Make /N=(ny,n) /O /R yval_set_inp
Make /N=(ny,ny) /O /R yval_ro
if(nx)
	Make /N=(nx) /O /T dname_set_inp
	Make /N=(nx,n) /O /R dval_set_inp
else
	Make /N=0 /O /T dname_set_inp
	Make /N=0 /O /R dval_set_inp
endif
//
if(ns<0)
	status_set_inp=1
endif
if(ns>=0)
	status_set_inp=str2num(inp0[p+1][ns+1])
endif
ycomp_set_inp=inp0[p+1][0]
Variable i=0,mx=0,my=0
do
	if(strsearch(kdes[i],".y",0)>0)		
		yname_set[my]=inp0[0][i+1]
		yval_set_inp[my][]=str2num(inp0[q+1][i+1])
		my+=1
	else
		if(strsearch(kdes[i],"status",0)<0)
			dname_set_inp[mx]=inp0[0][i+1]
			dval_set_inp[mx][]=str2num(inp0[q+1][i+1])
			mx+=1
		endif
	endif
i+=1
while(i<m)
KillWaves inp0,kdes
Variable /G table_inp=1
FindCorrespondence(opt)
avg_std()
end

Proc FindCorrespondence(opt)
Variable opt
Variable /G table_inp
Variable n=numpnts(comp_set)
duplicate /O comp_set comp_set_
duplicate /O smi_set,smi_set_
if(exists("class_set")==1)
	duplicate /O class_set class_set_
	duplicate /O full_set full_set_
	duplicate /O iupac_set iupac_set_
endif
Make /N=(n) /O /B status_set,status_set_
status_set=1
if(!table_inp)
	abort "no table input"
endif
Variable score=find_des_inp(opt)
if(score)
	printf "total of "+num2str(score)+" missing y's\r\n"
endif
duplicate /O dname_set_ dname_set
duplicate /O comp_set_ comp_set
duplicate /O smi_set_ smi_set
if(exists("class_set")==1)
	duplicate /O class_set_ class_set
	duplicate /O full_set_ full_set
	duplicate /O iupac_set_ iupac_set
	do_class_set()
endif
duplicate /O status_set_ status_set
duplicate /O dval_set_ dval_set
duplicate /O yval_set_ yval_set
duplicate /O dname_set_ dname_set
KillWaves yval_set_,dval_set_,dname_set_,comp_set_,smi_set_,status_set_
printf ">>> created new unshifted work set\r\n"
end

		
function find_des_inp(opt)
Variable opt
Wave /T dname_set,dname_set_
Wave /T comp_set_,comp_set
Wave /T smi_set,smi_set_
Wave /T class_set,class_set_
Wave /T full_set,full_set_
Wave /T iupac_set,iupac_set_
Wave /T ycomp_set_inp
Wave yval_set_inp,dval_set_inp
Wave /B status_set_inp
Wave /T yname_set,dname_set_inp
Variable n=numpnts(comp_set)
Variable m=numpnts(dname_set)
Variable /G ny=numpnts(yname_set)
Variable /G nx=numpnts(dname_set_inp)
Variable /G n_=numpnts(ycomp_set_inp)
Variable /G ndset=nx+m
Variable i,j,k
String s
//
Make /N=(ndset) /O /T dname_set_
for(k=0;k<m;k+=1)
	dname_set_[k]=dname_set[k]
endfor
//
if(nx)
	for(k=0;k<nx;k+=1)
		s=dname_set_inp[k]
		if(strsearch(s,"*",0)<0)
			s="*"+s
		endif
		dname_set_[m+k]=s
	endfor
endif
//
Make /N=(n) /O /B hash_comp
hash_comp=1
Variable i1=0
for(j=0;j<n_;j+=1)
	i=comp_match(j)
	if(i>=0)
		hash_comp[i]=0
		i1+=1
	endif
endfor
//
Variable score=sum(hash_comp,0,n-1)
if(opt)
	i1+=score
endif
Make /N=(ndset,i1) /O /R dval_set_
Make /N=(ny,i1) /O /R yval_set_
Make /N=(i1) /O /B status_set_
Make /N=(i1) /O /T smi_set_,comp_set_,class_set_,full_set_,iupac_set_
smi_set_=""
status_set_=0
yval_set_=Nan
dval_set_=Nan
//
i1=0
for(j=0;j<n_;j+=1)
	i=comp_match(j)
	if(i>=0)
		insert_des_inp(i1,i,j)
		i1+=1
	endif
endfor
//
if(opt)
	for(i=0;i<n;i+=1)
		if(hash_comp[i])
			insert_des_inp(i1,i,-1)
			i1+=1
		endif
	endfor
endif
return score
end


function insert_des_inp(i1,i,j)
Variable i1,i,j
Wave /T comp_set_,smi_set_,comp_set,smi_set,dname_set
Wave /T class_set, class_set_
Wave /T full_set, full_set_
Wave /T iupac_set, iupac_set_
Wave /B status_set_,status_set_inp
Wave yval_set_,yval_set_inp
Wave /T yname_set
Wave dval_set_,dval_set_inp,dval_set
Variable m=numpnts(dname_set)
Variable /G ny,nx
Variable k
//
comp_set_[i1]=comp_set[i]
smi_set_[i1]=smi_set[i]	
if(exists("class_set")==1)
	class_set_[i1]=class_set[i]
	full_set_[i1]=full_set[i]
	iupac_set_[i1]=iupac_set[i]
endif
for(k=0;k<m;k+=1)
	dval_set_[k][i1]=dval_set[k][i]
endfor
if(j>=0)
	 status_set_[i1]=status_set_inp[j]
	 for(k=0;k<ny;k+=1)
		yval_set_[k][i1]=yval_set_inp[k][j]
	endfor
	if(nx)
		for(k=0;k<nx;k+=1)
			dval_set_[m+k][i1]=dval_set_inp[k][j]
		endfor
	endif
endif
end

function comp_match(i)
Variable i
Wave /T comp_set,ycomp_set_inp
String s=ycomp_set_inp[i]
	return WhereInSWave(comp_set,s)
end

function avg_std()
Wave dval_set
Wave /T dname_set
Variable n=DimSize(dval_set,0)
Variable m=DimSize(dval_set,1)
Make /N=(n) /O /R av_set,std_set,pcvar_set
Variable i,j,x,q
av_set=0
std_set=0
for(i=0;i<n;i+=1)
	q=0
	for(j=0;j<m;j+=1)
		x=dval_set[i][j]
		if(!numtype(x))
			av_set[i]+=x
			std_set[i]+=x*x
			q+=1
		endif
	endfor
	av_set[i]/=q
	std_set[i]=sqrt(std_set[i]/q-av_set[i]^2)
endfor
pcvar_set=100*std_set/abs(av_set)
end



Window Table_variables_input() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,42.5,456,476) ycomp_set_inp,status_set_inp,dval_set_inp,yval_set_inp
	ModifyTable width(Point)=30,alignment(ycomp_set_inp)=0,rgb(ycomp_set_inp)=(0,0,65280)
	ModifyTable width(status_set_inp)=33,rgb(status_set_inp)=(65280,0,0),alignment(dval_set_inp)=0
	ModifyTable width(dval_set_inp)=71,alignment(yval_set_inp)=0,width(yval_set_inp)=53
	ModifyTable rgb(yval_set_inp)=(65280,0,0),elements(yval_set_inp)=(-3,-2)
EndMacro

Window Table_descriptor_set() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(578,46,928,375) dname_set,dval_set
	ModifyTable format(Point)=1,width(Point)=23,alignment(dname_set)=0,width(dname_set)=126
	ModifyTable rgb(dname_set)=(0,0,52224),alignment(dval_set)=0,format(dval_set)=3
	ModifyTable width(dval_set)=59
EndMacro

Window Table_des_stat() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(7,47,395,531) dname_set,av_set,std_set,pcvar_set
	ModifyTable format(Point)=1,width(Point)=0,alignment(dname_set)=0,width(dname_set)=108
	ModifyTable rgb(dname_set)=(1,4,52428),format(av_set)=3,digits(av_set)=4,sigDigits(av_set)=4
	ModifyTable format(std_set)=3,digits(std_set)=4,format(pcvar_set)=3,digits(pcvar_set)=1
EndMacro


Window Table_comp_set() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5,43,685,526) comp_set,smi_set,status_set,dval_set,yval_set
	ModifyTable format(Point)=1,width(Point)=21,size(comp_set)=9,alignment(comp_set)=0
	ModifyTable width(comp_set)=57,rgb(comp_set)=(0,0,65280),size(smi_set)=9,alignment(smi_set)=0
	ModifyTable width(smi_set)=314,width(status_set)=33,rgb(status_set)=(65280,0,0)
	ModifyTable alignment(dval_set)=0,format(dval_set)=3,width(dval_set)=62,elements(dval_set)=(-3,-2)
	ModifyTable alignment(yval_set)=0,format(yval_set)=3,width(yval_set)=62,rgb(yval_set)=(65280,0,0)
	ModifyTable elements(yval_set)=(-3,-2)
EndMacro


Window Table_variables_set() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(27.75,85.25,784.5,488.75) comp_set,smi_set,status_set,yval_set
	ModifyTable width(Point)=21,size(comp_set)=9,alignment(comp_set)=0,width(comp_set)=57
	ModifyTable rgb(comp_set)=(0,0,65280),size(smi_set)=9,alignment(smi_set)=0,width(smi_set)=314
	ModifyTable width(status_set)=32,rgb(status_set)=(65280,0,0),alignment(yval_set)=0
	ModifyTable format(yval_set)=3,width(yval_set)=62,elements(yval_set)=(-3,-2)
EndMacro


Window Table_variables_corr() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,42.5,231.75,185) yval_ro
	ModifyTable width(Point)=18,format(yval_ro)=3,width(yval_ro)=57
EndMacro

Window Table_comp_short() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(35.25,46.25,714.75,569.75) comp_set,smi_set,status_set,yval_set
	ModifyTable width(Point)=21,size(comp_set)=9,alignment(comp_set)=0,width(comp_set)=57
	ModifyTable rgb(comp_set)=(0,0,65280),size(smi_set)=9,alignment(smi_set)=0,width(smi_set)=314
	ModifyTable width(status_set)=33,rgb(status_set)=(65280,0,0),alignment(yval_set)=0
	ModifyTable format(yval_set)=3,width(yval_set)=62,rgb(yval_set)=(0,0,65280),elements(yval_set)=(-3,-2)
EndMacro


Window Table_descriptors_short() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(474.75,220.25,636.75,473.75) dname_set
	ModifyTable width(Point)=24,alignment(dname_set)=0,width(dname_set)=105,rgb(dname_set)=(0,0,65280)
EndMacro


Window Table_variables_short() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(796,381,1002,913) yname_set,r2y
	ModifyTable format(Point)=1,width(Point)=30,alignment(yname_set)=0,width(yname_set)=94
	ModifyTable rgb(yname_set)=(0,0,52224),format(r2y)=3,digits(r2y)=4,width(r2y)=58
EndMacro

//====================================================================

Proc HashSelection()
GetSelection table,$"Table_descriptors_short",1
dname_set[V_StartRow,V_EndRow]=hash_(dname_set[p])
end

Proc UnhashSelection()
GetSelection table,$"Table_descriptors_short",1
dname_set[V_StartRow,V_EndRow]=replacestring("#",dname_set[p],"")
end

Proc StarSelection()
GetSelection table,$"Table_descriptors_short",1
dname_set[V_StartRow,V_EndRow]=star_(dname_set[p])
end

Proc HashSmallStd(threshold)
Variable threshold=5
Prompt threshold,"threshold % of average"
PauseUpdate
Silent 1
avg_std()
String s
Variable i=0,n=numpnts(pcvar_set)
do
	s=dname_set[i]
	if ((pcvar_set[i]<threshold)) //&& (strsearch(s,"#",0)<0))
		dname_set[i] = hash_(s)
	endif
i+=1
while(i<n)
end


Proc HashAll3D()
hash_G98_descriptors()
dname_set=hash_3D(dname_set[p])
end

Proc HashAll2n3D()
hash_G98_descriptors()
dname_set=hash_3D(dname_set[p])
dname_set=hash_2D(dname_set[p])
end

Proc UnstarAllDescriptors()
	dname_set=replacesymbols(dname_set,"*","")
end

Proc UnstarAllYs()
	yname_set=replacesymbols(yname_set,"*","")
end

Proc StarAllYs()
	yname_set=addsymbolcheck(yname_set,"*",0)
end

function /S addsymbolcheck(s,symb,opt)
String s,symb
Variable opt
s=replacesymbols(s,symb,"")
	if(opt)
		return s+symb
	else
		return symb+s
	endif
end

Proc UnhashAllDescriptors()
	dname_set=replacesymbols(dname_set,"#","")
end


Macro IdentifyCompound()
Variable a=pcsr(A)
printf ">>> %s :: %s\r\n", comp_set[a], smi_set[a]
if(exists("class_set")==1)
	printf ">>> %s : %s : %s", class_set[a], full_set[a], iupac_set[a]
endif
end

function hash_G98_descriptors()
Wave /T dname_set
Variable i,n=numpnts(dname_set)
String s
do
	s=hash_(dname_set[i])
	if (strsearch(s,".G98",0)>=0) 
		dname_set[i] = hash_(s)
	endif
i+=1
while(i<n)
end

function /S hash_(s)
String  s
//s=replacesymbols(s,"\r","")
//s=replacesymbols(s,"\n","")
s=replacesymbols(s," ","")
if(strsearch(s,"#",0)>=0)
	return s
else
	return "#"+s
endif
end

function /S star_(s)
String  s
//s=replacesymbols(s,"\r","")
//s=replacesymbols(s,"\n","")
s=replacesymbols(s," ","")
if(strsearch(s,"*",0)>=0)
	return s
else
	return "*"+s
endif
end


function /S hash_3D(s)
String  s
//s=replacesymbols(s,"\r","")
//s=replacesymbols(s,"\n","")
s=replacesymbols(s," ","")
if(strsearch(s,"#",0)>=0)
	return s
endif
if(strsearch(s,"Geom",0)>=0)
	return "#"+s
endif
Wave /T dname_m,type_m
Variable m=WhereInSWave(dname_m,s)
if((m>=0) && (strsearch(type_m[m],"3D",0)>=0))
	return "#"+s
endif
return s
end

function /S hash_2D(s)
String  s
String q,list1D="n;Num;Ns;Na;Nd;Ss;Sa;Sd;MAX;MIN"
//s=replacesymbols(s,"\r","")
//s=replacesymbols(s,"\n","")
s=replacesymbols(s," ","")
if(strsearch(s,"#",0)>=0)
	return s
endif
if(strsearch(s,"Geom",0)>=0)
	return "#"+s
endif
Wave /T dname_m,type_m
Variable m=WhereInSWave(dname_m,s)
if((m>=0) && (strsearch(type_m[m],"2D",0)>=0))
	return "#"+s
endif
Variable i,n=ItemsInList(list1D),flag=0
for(i=0;i<n;i+=1)
	q=StringFromList(i,list1D)
	if(!strsearch(s,q,0))
		flag=1
	endif
endfor
if(!flag)
  s="#"+s
endif
return s
end


// creating sets by descriptor import ==========================================

Proc CreateSetByList_des(suffix)
String suffix="_out"
Prompt suffix,"suffix: will try it first"
Silent 1
PauseUpdate
String /G y_transform=""
Variable /G ny=0
Make /N=0 /R /O yval_set
Make /N=0 /T /O yname_set
String /G s_filename,s_path
Variable /G nset,ndset
String s
LoadWave/Q/J/N=wave/O/K=2
String /G s_path_data=s_path
String /G s_path_des=s_path+"des:"
newpath /O data s_path_data
newpath /O des s_path_des
String /G lib=removefileextension(s_filename)
duplicate /O wave0 comp_set
duplicate /O wave1 smi_set
KIllWaves wave0,wave1
Variable n=numpnts(comp_set)
Make /N=(0,n) /O /R dval_set
String /G directory_list=IndexedFile(des,-1,".des")
test_des_set_list(suffix)
printf "%d compounds x %d descriptors\r\n",nset,ndset
duplicate /O comp_ comp_set
duplicate /O smi_ smi_set
KillWaves comp_,smi_
Make /N=(ndset,nset) /O /R dval_set
Make  /N=(ndset) /O /T dname_set
make_set_des() 
ExportSet(lib)
Variable /G table_inp=0
end

Proc CreateSetAll_des(library)
String library="set"
Prompt library,"library name"
Silent 1
PauseUpdate
String /G y_transform=""
Variable /G ny=0
Make /N=0 /R /O yval_set
Make /N=0 /T /O yname_set
String /G s_filename,s_path
Variable /G nset,ndset
String s
Open /R /T=".des" refnum 
Close refnum
removefilepath_(s_filename,":")
String /G s_path_des=s_path
newpath /O des s_path_des
String /G lib=library
String /G directory_list=IndexedFile(des,-1,".des")
test_des_set_all()
printf "%d compounds x %d descriptors\r\n",nset,ndset
duplicate /O comp_ comp_set
duplicate /O smi_ smi_set
KillWaves comp_,smi_
Make /N=(ndset,nset) /O /R dval_set
Make  /N=(ndset) /O /T dname_set
make_set_des() 
ExportSet(lib)
Variable /G table_inp=0
end

Proc Add_Classes2Descriptors()
add_class_des()
end

function add_class_des()
String /G class_list
Wave /T class_set, dname_set
Wave dval_set
Variable i,j
String s
Variable n=numpnts(dname_set)
Variable m=numpnts(class_set)
Variable k=ItemsInList(class_list)
//
Make /N=(n+k) /O /T dname_set_
dname_set_= dname_set
for (j=0;j<k;j+=1)
	s=StringFromList(j,class_list)
	dname_set_[n+j]="in_"+s+".par"
endfor
//
Make /N=(n+k,m) /O /R dval_set_
dval_set_=0
dval_set_[0,n-1][] = dval_set[p][q]
for(i=0;i<m;i+=1) 
	j=WhichListItem(class_set[i],class_list)
	dval_set_[n+j][i] = 1
endfor
duplicate /O /T dname_set_ dname_set
duplicate /O dval_set_ dval_set
KillWaves /Z dname_set_, dval_set_
avg_std()
end

function do_class_set()
Wave /T class_set
Wave status_set
Variable i, j, m=1, n=numpnts(class_set)
String /G class_list=""
Make /O /N=(n) /R color_set
String s
for(i=0;i<n;i+=1)
	if(status_set[i]<2)
		s = TrimString(class_set[i])
		j = WhichListItem(s, class_list)
		if(j<0)
			class_list+=s+";"
		endif
	endif
endfor
class_list = SortList(class_list,";",16)
color_set =  WhichListItem(class_set, class_list)
printf "\r%d classes :: %s",m,class_list
color_index(StringFromList(0,class_list))
// histograms
m=ItemsInList(class_list)
Make /N=(m) /T /O classes
classes= StringFromList(p, class_list)
Make /N=(m) /R /O class_total, class_tested
class_total=0
class_tested=0
for(i=0;i<n;i+=1)
	s = TrimString(class_set[i])
	j = WhichListItem(s, class_list)
	if(j<0)
		continue
	endif
	class_total[j]+=1
	if(status_set[i]==1)
		class_tested[j]+=1
	endif
endfor
class_total*=100/n
class_tested*=100/n
end

function color_index(s)
String s
Wave /T class_set
Variable i, n=numpnts(class_set)
Make /O /N=(n) /R index_set
index_set=0
for(i=0;i<n;i+=1)
	if(stringmatch(s,class_set[i]))
		index_set[i]=1
	endif
endfor
end


Proc CreateSet_HDF5_des(keep)
Variable keep=0
Prompt keep,"0 - new, 1 - keep adding"
Silent 1
Variable /G refnum
HDF5OpenFile /R /I refnum as ""
read_des_h5(keep)
avg_std()
HDF5CloseFile refnum
String /G lib=""
Variable /G table_inp=0
end

function read_des_h5(keep)
Variable keep
Wave /T dname_set
String s,q,list
Variable u,i,i1,j,m,n,good=0,bad=0
Variable /G refnum,ndset,nset,ny
String /G file_set,smiles_set,id_set
Make /N=0 /T /O sdummy,dname_type
//
if(!keep)
	id_set=""
	file_set=""
	smiles_set=""
	nset=0
endif
//
HDF5ListGroup /F refnum, "/"
list=replacesymbols(S_HDF5ListGroup,"/names;","")
n=ItemsInList(list)
HDF5LoadData /Q/O /N=dname_type refnum,"/names"
Variable nd=numpnts(dname_type)
Make /N=(nd) /O /R dval,av_set,std_set
Make /N=(nd) /O /B hash_set
av_set=0; std_set=0; hash_set=0
//
for(i=0;i<n;i+=1)
	s=StringFromList(i,list)
	HDF5LoadData /Q/O /N=dval refnum, s
	HDF5LoadData /Q /O /A="smiles" /N=sdummy refnum, s
	s=s[1,strlen(s)-1]		
	q=sdummy[0]
	u=WhichListItem(q, smiles_set)
	if(u<0)
		id_set+=s+";"
		smiles_set+=q+";"
		good+=1
	   av_set+=dval
	   std_set+=dval^2	
	else
		if(strlen(smiles_set))
			printf "\r\n(%3d) %s : %s : same as (%3d) %s",bad,s,q,u,StringFromList(u,id_set)
			bad+=1
		endif
	endif
endfor
//
if(!keep)
	make_log_des(good)
else
	Make /N=(ndset) /O dname_pos
	dname_pos=WhereInSWave(dname_type,dname_set[p])
endif
//
printf "\r\nHDF5: found %d compounds with %d descriptors",n,nd
if(!good)
	printf "\r\n>> no new smiles\r\n"
	return 0
else
	printf "\r\n>> reduced set to %d unique compounds with %d descriptors\r\n",good,ndset
	if(nset)
		printf "\r\n>> added to %d existing compounds\r\n",nset
	endif
endif
//
Make /N=(ndset,nset+good) /O /R dval_set
Make /N=(nset+good) /O /T comp_set,smi_set,class_set,full_set,iupac_set

for(i=0;i<good;i+=1)
	i1=nset+i
	s=StringFromList(i1,id_set)
	q="/"+s
	HDF5LoadData /Q/O /N=dval refnum,q
	
	sdummy=""
	HDF5LoadData /Q /O /A="smiles" /N=sdummy refnum,q
	if(WhichListItem(s,file_set)>=0)
		s="N"+s
	endif
	file_set+=s+";"
	comp_set[i1]=s		
	smi_set[i1]=sdummy[0]
	
	sdummy=""
	HDF5LoadData /Z /Q /O /A="class" /N=sdummy refnum,q
	class_set[i1]=sdummy[0]
	
	sdummy=""
	HDF5LoadData /Z /Q /O /A="full_name" /N=sdummy refnum,q
	full_set[i1]=sdummy[0]
	
	sdummy=""
	HDF5LoadData /Z /Q /O /A="IUPAC_name" /N=sdummy refnum,q
	iupac_set[i1]=sdummy[0]
	
	
	if(!keep)
		m=0	
		for(j=0;j<nd;j+=1)
			if(!hash_set[j])
				dval_set[m][i1]=dval[j]
				m+=1
			endif
		endfor
	else
		for(j=0;j<ndset;j+=1)
			m=dname_pos[j]
			if(m<0)
				dval_set[j][i1]=nan
			else
				dval_set[j][i1]=dval[m]
			endif
		endfor
	endif
endfor
//
if(!keep)
	Make  /N=(ndset) /O /T dname_set
	m=0
	for(j=0;j<nd;j+=1)
		if(!hash_set[j])
			dname_set[m]=dname_type[j]
			m+=1
		endif
	endfor
endif
//
nset+=good
Make /N=(nset) /O /B status_set
Make /N=(ny,nset) /O /R yval_set
yval_set=nan
status_set=0
KillWaves /Z sdummy,dname_pos
return good
end


Proc ParseSet_des(suffix,opt)
Variable opt=1
prompt opt,"if 1 - do IQR test, 2 - only IQR test"
String suffix=lib
Prompt suffix,"if empy,does not save and export"
Silent 1
PauseUpdate
String /G lib
Silent 1
PauseUpdate
simple_stat_test(opt)
avg_std()
if(strlen(suffix))
	lib=suffix
	SaveSet(suffix)
	ExportSet(lib+"_")
endif
end

Proc SaveSet(suffix)
String suffix=removefileextension(lib)
Silent 1
PauseUpdate
String /G lib=suffix
String s
Variable /G ny=numpnts(yname_set)
if(strlen(suffix))
	printf ">>> saving unshifted sets\r\n"

	s="dval_set;dname_set;comp_set;smi_set;av_set;std_set"
	CopyWaveList(s,suffix,1)
	if(ny)
		s="yval_set;yname_set;status_set;yav_set;ystd_set"
		CopyWaveList(s,suffix,1)
	endif
	mkucomp()
	mkparam()
endif
end


proc ExportSet(s)
String s=removefileextension(lib)
Silent 1
PauseUpdate
String /G lib=s
if(strsearch(s,"lib_",0)<0)
	s="lib_"+lib
endif
export_set_csv(s)
end


function Load_descriptors()
Variable opt
Variable /G refnum
Wave /T data0,dname
Wave dval
String s,q
Variable i,m=0,n
//
if(FindInFile("Descriptors:")==0)
	return 0
endif
FeedUntil("END")
n=numpnts(data0)
//
for(i=0;i<n;i+=1)
	s=data0[i]
	if(strsearch(s,"=",0)>=0)
		Make /N=(m+1) /O /T dname
		Make /N=(m+1) /O /R dval
		q=StringFromList(0,s,"=")
		dname[m]=replacesymbols(q,","," ")
		dval[m]=str2num(StringFromList(1,s,"="))
		m+=1
	endif
endfor	
//
String /G smiles=""
if(FindInFile("SMILE"))
	FeedUntil("END")
	n=numpnts(data0)
	for(i=0;i<n;i+=1)
		smiles+=data0[i]
	endfor
endif
end


function test_des_set_list(suffix)
String suffix
Wave /T comp_set,smi_set
String /G directory_list
Variable /G refnum
Variable i,n=numpnts(comp_set)
String /G file_set=""
Variable m
String s,q
Wave /T dname
Wave dval
Variable nd,k,k1,k2
String /G smiles
Variable /G nset=0,ndset=0
Variable /G reflog
Open /P=home reflog "log_des_test.csv"
for(i=0;i<n;i+=1)
	s=comp_set[i]
	if(strlen(s))
		m=FindInListExact(s+suffix+".des",directory_list,";")
		if(m<0)
			m=FindInListExact(s+".des",directory_list,";")
			if(m<0)
				continue
			endif
		endif
	else
		continue
	endif
	s=StringFromList(m,directory_list)
	print i,s
	Open/R /P=des refNum as s
	Load_descriptors()
	Close refnum
	if(!i)
		print ">>> created prototype"
		duplicate /O dname dname_type
		nd=numpnts(dname_type)
		Make /N=(nd) /O /R av_set,std_set
		Make /N=(nd) /O /B hash_set
		av_set=0; std_set=0; hash_set=0
	endif
	m=CompareSWaves(dname_type,dname)
	if(m)
		printf ">>> does not match the prototype\r\n"
		fprintf reflog,"mismatch for %s\r\n",s
		des2prototype()
	endif
	av_set+=dval
	std_set+=dval^2	
	Make /N=(nset+1) /O /T comp_,smi_
	comp_[nset]=comp_set[i]
	smi_[nset]=smi_set[i]
	file_set+=s+";"
nset+=1	
endfor
make_log_des(nset)
end


function test_des_set_all()
String /G directory_list
Variable /G refnum
Variable i,n=numpnts(comp_set)
String /G file_set=""
Variable m
String s,q
Wave /T dname
Wave dval
Variable nd,k,k1,k2
String /G smiles
Variable /G nset=0
Variable /G reflog
n=ItemsInList(directory_list)
if(!n)
	abort "no descriptors in the directory"
endif
Make /N=(n) /T /O comp_set,smi_set
for(i=0;i<n;i+=1)
	s=StringFromList(i,directory_list)
	print s
	comp_set[i]=replacesymbols(removefileextension(s),"_out","")	
	Open/R /P=des refNum as s
	Load_descriptors()
	smi_set[i]=smiles
	Close refnum
	if(!i)
		duplicate /O dname dname_type
		nd=numpnts(dname_type)
		Make /N=(nd) /O /R av_set,std_set
		Make /N=(nd) /O /B hash_set
		av_set=0; std_set=0; hash_set=0
	endif
	m=CompareSWaves(dname_type,dname)
	if(m)
		printf ">>> does not match the prototype\r\n"
		fprintf reflog,"mismatch for %s\r\n",s
		des2prototype()
	endif
	av_set+=dval
	std_set+=dval^2	
	Make /N=(nset+1) /O /T comp_,smi_
	comp_[nset]=comp_set[i]
	smi_[nset]=smi_set[i]
	file_set+=s+";"
nset+=1	
endfor
make_log_des(nset)
end


function make_log_des(n)
Variable n
Wave /B hash_set
Wave std_set,av_set
Wave /T dname_type
Variable /G ndset=0
Variable nd=numpnts(dname_type)
Variable m,i
String s,meaning
av_set/=n
std_set=sqrt(std_set/n-av_set*av_set)
variable /G reflog
Open /P=home reflog "log_des_test.csv"
fprintf reflog,"\r\n\r\nstatistics for input:\r\n\r\n#,descriptor,avg,st.dev,hash,meaning\r\n"
hash_set=0
for(i=0;i<nd;i+=1)		
	m=std_set[i]
	if( (m<1e-4) || (numtype(m)) )
		hash_set[i]=1		
	else
		ndset+=1
	endif
	s=dname_type[i]
	s=replacesymbols(s,","," ")
	meaning=des_explain_(s)
	fprintf reflog,"%d,%s,%.4f,%.4f,%d,,%s\r\n",i,s,av_set[i],std_set[i],hash_set[i],meaning
endfor
Close(reflog)
end

function des2prototype()
Wave /T dname,dname_type
Wave dval
Variable k,j,m=numpnts(dname_type)
Variable n=numpnts(dname)
String s=""
for(j=0;j<n;j+=1)
	s+=dname[j]+";"
endfor
duplicate /O dval dval_
Make /N=(m) /O /R dval
for(j=0;j<m;j+=1)
	k=WhichListItem(dname_type[j],s)
	//print j,dname_type[j],k
	if(k>=0)
		dval[j]=dval_[k]
	else
		dval[j]=nan
	endif
endfor
KillWaves dval_
Make /N=(m) /O /T dname
dname=dname_type
end

function export_set_csv(c)
String c
Wave /T dname_set,comp_set,smi_set,yname_set
Wave dval_set,yval_set
Wave /B status_set
Variable i,j,z
Variable ny=DimSize(yval_set,0)
Variable nd=DimSize(dval_set,0)
Variable n=DimSize(dval_set,1)
Variable /G refnum
Open /T="TEXT" /P=home refnum as c+".csv"
fprintf refnum,"%s,smiles",c
for(j=0;j<nd;j+=1)
	fprintf refnum,",%s",dname_set[j]
endfor
if(ny)
	for(j=0;j<ny;j+=1)
		fprintf refnum,",%s",yname_set[j]
	endfor
	fprintf refnum,",status"
endif
//
fprintf refnum,"\r\n"
for(i=0;i<n;i+=1)
	fprintf refnum,"%s,%s",comp_set[i],smi_set[i]
	for(j=0;j<nd;j+=1)
		fprintf refnum,",%.6f",dval_set[j][i]
	endfor
	if(ny)
		for(j=0;j<ny;j+=1)
			z=yval_set[j][i]
			if(numtype(z))
				fprintf refnum,",Nan"
				continue
			endif	
			if(!z)
				fprintf refnum,",0"
				continue
			endif
			if(abs(z)>1e-2)
				fprintf refnum,",%.6f",z
				continue
			endif
			fprintf refnum,",%.6e",z
		endfor
		fprintf refnum,",%d",status_set[i]
	endif
	fprintf refnum,"\r\n"
endfor
Close refnum
end

function make_set_des()
Variable opt
Variable /G refnum
String /G file_set
Wave /T comp_set,dname,smi_set,dname_set,dname_type
Wave dval,dval_set
Wave /B hash_set
Variable z,i,j,m
String s,fs
String /G directory_list
Variable n=numpnts(comp_set)
Variable nd=numpnts(dname_type)
for(i=0;i<n;i+=1)
	fs=StringFromList(i,file_set)
	Open/R /P=des refNum as fs
	Load_descriptors()
	if(CompareSWaves(dname_type,dname))
		des2prototype()
	endif
	Close refnum
	printf ">>> %s\r\n",fs	
	m=0	
	for(j=0;j<nd;j+=1)
		if(!hash_set[j])
			dval_set[m][i]=dval[j]
			m+=1
		endif
	endfor
endfor
m=0
for(j=0;j<nd;j+=1)
	if(!hash_set[j])
		dname_set[m]=dname_type[j]
		m+=1
	endif
endfor
end

function value_set_des(dn,dv,s,i)
Wave /T dn
Wave dv
String s
Variable i
if(i>=DimSize(dv,1))
	return nan
endif
Variable j=WhereInSWave(dn,s)
if(j>0)
	return dv[j][i]
else
	return nan
endif
end


//========= combines "anion" and "cation" sets

Proc Export_salts()
LoadWave/Q/J/N=wave/O/K=2
duplicate /O wave0 cation_list
duplicate /O wave1 anion_list
KillWaves wave0,wave1
export_salt_set()
end

Proc Test_salts()
Silent 1
PauseUpdate
LoadWave/Q/J/N=wave/O/K=2/V={","," $",0,1}
duplicate /O wave0 cation_list
duplicate /O wave1 anion_list
KillWaves wave0,wave1
Variable i=0,n=numpnts(cation_list),flag
Variable /G refnum
Variable m=0
do
	flag=0
	if(WhereInSWave(comp_set_cation,cation_list[i])<0)
		printf "line %d: cannot find cation %s in the cation set\r\n",i+1,cation_list[i]
		flag=1
	endif
	if(WhereInSWave(comp_set_anion,anion_list[i])<0)
		printf "line %d: cannot find anion %s in the anion set\r\n",i+1,anion_list[i]
		flag=1
	endif
	if(flag)
		m+=1
	endif
i+=1
while(i<n)
if(m)
	printf ">>> %d salts not in the cation or anion sets\r\n",m
else
	printf ">>> all salts are in the cation and anion sets\r\n"
endif
end


function export_salt_set()
Variable /G refnum
Wave /T cation_list,anion_list
Wave /T dname_set_anion
Wave /T dname_set_cation
Wave /T smi_set_anion
Wave /T smi_set_cation
Wave /T comp_set_anion
Wave /T comp_set_cation
Wave  dval_set_anion
Wave  dval_set_cation
Variable i,j,ia,ic
String c,a,s
Open /T="TEXT" /P=home refnum "lib_salts.csv"
Variable n=numpnts(cation_list)
Variable m=0
Variable na=DimSize(dval_set_anion,0)
Variable nc=DimSize(dval_set_cation,0)
fprintf refnum,"salts,smiles (+ -)"
//
// salt specific descriptors
//
fprintf refnum,",Eac.salt"
//
for(j=0;j<nc;j+=1)
	fprintf refnum,",%s(+)",dname_set_cation[j]
endfor
for(j=0;j<na;j+=1)
	fprintf refnum,",%s(-)",dname_set_anion[j]
endfor
fprintf refnum,"\r\n"
//
Variable ra,rc,Eac
for(i=0;i<n;i+=1)
	c=cation_list[i]
	a=anion_list[i]
	sprintf s,"[%s+][%s-]",c,a
	printf "%3d:\t%s\r\n",i+1,s
	if( !(strlen(a)+strlen(c)) )
		printf ">>> missing input\r\n"
		continue
	endif
	ic=WhereInSWave(comp_set_cation,c)
	if (ic<0) 
		printf ">>> missing cation\n\r"
		continue
	endif
	ia=WhereInSWave(comp_set_anion,a)
	if (ia<0) 
		printf ">>> missing anion\n\r"
		continue
	endif
	fprintf refnum,"%s,%s   %s",s,smi_set_cation[ic],smi_set_anion[ia]
	//
	// salt specific descriptors
	//
	rc=value_set_des(dname_set_cation,dval_set_cation,"*a0_SCRF.G98",ic)
	ra=value_set_des(dname_set_anion,dval_set_anion,"*a0_SCRF.G98",ia)
	Eac=100/(ra+rc)
	fprintf refnum,",%.3f",Eac
	//		
	for(j=0;j<nc;j+=1)
		fprintf refnum,",%.3f",dval_set_cation[j][ic]
	endfor
	for(j=0;j<na;j+=1)
		fprintf refnum,",%.3f",dval_set_anion[j][ia]
	endfor
	fprintf refnum,"\r\n"
m+=1
endfor
Close refnum
printf "\r\nexported %d salts out of %d in the list\r\n",m,n
end


//====================== stat tests ===================


Proc simple_stat_test(opt)
Variable opt
stat_test_dval_set(opt,0.9)
reduce_by_hash(1)
duplicate /O dval_set_ dval_set
duplicate /O dname_set_ dname_set
KillWaves dval_set_,dname_set_,hash_set
end

Proc std_stat_test()
std_test_dval_set()
Variable /G nd_reduced=reduce_by_hash(0)
if(nd_reduced)
	duplicate /O dval_set_ dval_set
	duplicate /O dname_set_ dname_set
	KillWaves dval_set_,dname_set_
endif
KillWaves hash_set
end

function std_test_dval_set()
Wave dval_set
Wave std_set
Wave /T dname_set
Variable nd=DimSize(dval_set,0)
Make /N=(nd) /B /O hash_set
hash_set=0
Variable i,j,m
m=0
for(i=0;i<nd;i+=1)
	if(std_set[i]<1e-6)
			hash_set[i]=1
			m+=1
			printf "(%d) std fail for %s\r\n",i+1,dname_set[i]
		endif
	endfor
if(m)
	printf "std test: rejected %d out of %d descriptors\r\n",m,nd
endif
end

function stat_test_dval_set(opt,ro0)
Variable ro0,opt
Wave dval_set
Wave /T dname_set
Variable nd=DimSize(dval_set,0)
Make /N=(nd) /B /O hash_set
hash_set=0
Variable i,j,ro,m
// IQR test, options 1 and 2
if(opt)
	m=0
	for(i=0;i<nd;i+=1)
		if(IQR_2D(dval_set,i,1)==0)
			hash_set[i]=2
			m+=1
			printf "(%d) IQR fail for %s\r\n",i+1,dname_set[i]
		endif
	endfor
	printf "IQR test: rejected %d out of %d descriptors\r\n",m,nd
endif
// correlation test, options 0 and 2
if(opt!=2)
	m=0
	for(i=0;i<nd;i+=1)
		if(hash_set[i])
			continue
		endif
		for(j=i+1;j<nd;j+=1)
			if(hash_set[j])
				continue
			endif
			ro=Correlation_2D(dval_set,i,j,1)
			if(ro>ro0)
				hash_set[j]=3
				m+=1
				printf "(%d,%d) ro=%.3f :: %s  vs.  %s\r\n",i+1,j+1,ro,dname_set[i],dname_set[j]
			endif
		endfor
	endfor
	printf "correlation test (ro>%.2f): rejected %d ot of %d descriptors\r\n",ro0,m,nd
endif
end

function reduce_by_hash(opt) // 1 - ignore starred
Variable opt
Wave dval_set
Wave /B hash_set
Wave /T dname_set
Variable nd=DimSize(dval_set,0)
Variable n=DimSize(dval_set,1)
Variable i,j,m
m=0
for(j=0;j<nd;j+=1)
	if(!hash_set[j])
		m+=1
	else // make exceptions for * descriptors
		if((opt) && (keep_des(j,dname_set,"*")))
			hash_set[j]=0
			m+=1
		endif
	endif
endfor
if(m!=nd)
	printf "accepted %d out of %d descriptors\r\n",m,nd
endif
//
	Make /N=(m,n) /O /R dval_set_
	Make /N=(m) /T /O dname_set_
	m=0
	for(j=0;j<nd;j+=1)
		if(!hash_set[j])
			dname_set_[m]=dname_set[j]
			m+=1
	endif
	endfor
	for(i=0;i<n;i+=1)
		m=0
		for(j=0;j<nd;j+=1)
			if(!hash_set[j])
				dval_set_[m][i]=dval_set[j][i]
				m+=1
			endif
		endfor
	endfor
return m
end

function keep_des(j,w,symb)
Variable j
String symb
Wave /T w
Variable i,n=ItemsInList(symb)
String q,s=w[j]
for(i=0;i<n;i+=1)
	q=StringFromList(i,symb)
	if(strsearch(s,q,0)>=0)
		return 1
	endif
endfor
return 0
end




//======================= MM linear regression analyses ===============
//
// genetic algorithm reduction
//

Proc Prepare_MRL_genetics()
Silent 1
PauseUpdate
Variable /G show_comp_index=0
String /G show_comp=""
avg_std()
Variable /G Ng=min(Ng,floor(DimSize(dval_set,0)*.6))
printf "%d random descriptors\r\n",Ng
Make /N=0 /O /R yp1,yv1,ep1
Make /O /N=2 line45
Make /N=0 /O /R xv,yv,xvb,yvb,yav,ystd,yp,ypb
Make /N=0 /O M_product,M_A,M_B
mkucomp()
mkparam()
choose_y(yname_set)
Variable /G yshow=CheckCircularBelong(ysubset,yshow)
String   /G yname_show=yname_set[yshow]
//
Prepare_regression()
Prepare_Genetics(5e4)
DoWindow /K Table_coeff
DoWindow /K GeneticsPanel
GeneticsPanel()
end

Proc Prepare_regression()
Variable /G r2show=nan
Variable /G solution_found=0
	subset_train_valid()
	choose_fixed_subset()
	duplicate /O status_set color_set
end


Proc Prepare_Genetics(NTmax)
Variable NTmax
//
	Variable /G DBG=0
	Variable /G Nv=50
	Variable /G Nvs=10
	Variable /G Nst=2
	Variable /G NT=NTmax
	Variable /G rmsB=0
	Variable /G IT=0
//
end

Macro MLR_Genetics(n,l)
Variable n=Ng,l=0
Prompt n,"number of genes"
Prompt l,"regularization, lambda (=0 by default)"
Silent 1
PauseUpdate
Variable /G do_regression=1
Variable /G lambda=l,lev=nan
Variable /G reg_type=0  // linear regression
Variable /G Ng=n
Variable /G counter_label = 0
Prepare_MRL_genetics()
make_random_vector()
Genetic(NT)
hat_eigen_MLR()
end

Proc CheckFixedDes()
choose_fixed_subset()
end

function choose_y(yname_set)
Wave /T yname_set
Variable n=numpnts(yname_set)
Make /N=(n)  /O /B yhash_set
yhash_set=keep_des(p,yname_set,"*")
Variable m=sum(yhash_set,0,n-1)
Make /N=(m) /O /I ysubset
Variable i
m=0
for(i=0;i<n;i+=1)
	if(yhash_set[i])
		ysubset[m]=i
		m+=1
		printf ">>> %s: keep\r\n",yname_set[i]
	endif
endfor
end

function Check_Nan_y(n)
Variable n
Wave /I ysubset
Wave yval_set
Variable i,m=numpnts(ysubset)
Variable k,u
	for(i=0;i<m;i+=1)
		k=ysubset[i]
		u=yval_set[k][n]
		if(numtype(u))
			return 0
		endif
	endfor
return 1
end

function subset_train_valid()
Wave /B status_set
Variable i,kt,kv,n=numpnts(status_set)
kt=0;kv=0
for(i=0;i<n;i+=1)
	if( Check_Nan_y(i) )
		if (status_set[i]==1) 
			Make /N=(kt+1) /O /I subset_train
			subset_train[kt]=i
			kt+=1
			continue
		endif
	else
		if(status_set[i]==1)	
			status_set[i]=2
		endif
	endif
	if(status_set[i])	
		Make /N=(kv+1) /O /I subset_valid
		subset_valid[kv]=i
		kv+=1
	endif
endfor
end

function make_random_vector()
Variable /G Ng
Wave /I subset
Make /N=(Ng) /O /R VB
choose_random_subset()
VB=subset
make_V_comment(VB)
end

function make_V_comment(v)
Wave v
Variable k=numpnts(v)
Wave /T dname_set
Make /N=(k) /O /T V_comment
V_comment=dname_set[v[p]]
end

function make_VB_comment()
String /G nindex
Variable k=ItemsInList(nindex)
Make /N=(k) /O /T VB_comment
VB_comment=index2name(StringFromList(p,nindex))
end

function choose_fixed_subset()
Wave /T dname_set
Variable m=DimSize(dname_set,0)
Make /N=(m) /O /B hash_set,sub_set
hash_set=keep_des(p,dname_set,"*") // * are always in the set"
Variable n=sum(hash_set,0,m-1)
Make /N=(n) /O /I subset_fixed
Variable j
n=0
for(j=0;j<m;j+=1)
	if(hash_set[j])
		printf "(%d) %s: keep\r\n",j,dname_set[j]
		subset_fixed[n]=j
		n+=1
	endif
endfor
if(n)
	printf "%d fixed descriptors\r\n",n
endif
Variable /G nfixed=n
//check_std_set()
hash_set+=keep_des(p,dname_set,"#") // # are always excluded
return n
end

function choose_random_subset()
Variable /G Ng
Wave hash_set,sub_set
Make /N=(Ng) /O /I subset
sub_set=hash_set
random_mask(sub_set,subset,Ng,0)
end



function random_mask(mask,list,k,n)
Wave /B mask
Wave /I list
Variable k,n
Variable i,u,z,j
//
if(!k)
	return 0
endif
//
Variable m1=0,m=numpnts(mask)
for(i=0;i<m;i+=1)
	if(mask[i]==0)
		m1+=1
	endif
endfor
//
for(i=0;i<k;i+=1)
	u=irand(m1)
	z=0
	for(j=0;j<m;j+=1)
		if(mask[j]==0)
			if(z==u)
				mask[j]=1
				list[n]=j
				n+=1
				m1-=1
				break
			endif
			z+=1
		endif
	endfor
endfor
Sort list,list
return n
end


function point_mutation(v,nm)
Wave v
Variable nm
Variable /G Ng
Wave hash_set,sub_set
Wave /I subset
Variable i,n=0,j,z=nm/Ng
sub_set=hash_set
sort v,v
for(i=0;i<Ng;i+=1)
	if(rand()>z)
		j=v[i]
		sub_set[j]=1
		subset[n]=j
		n+=1
	endif
endfor
random_mask(sub_set,subset,Ng-n,n)
v=subset
return Ng-n
end

function EvVect(v)
Wave v
	return  EvVect_(v,0)
end

function EvVect_(v,opt)
Wave v
Variable opt
Wave subset_fixed,dval_set
Wave xv,yv
Variable k=numpnts(v)
Variable /G nfixed
Variable m=nfixed+k
Make /N=(m) /O /R VB_
if(nfixed)
	VB_[0,nfixed-1]=subset_fixed[p]
endif
VB_[nfixed,m-1]=v[p-nfixed]
//
Variable /G reg_type
if(reg_type==0) 
	return EvSubset_MLR(VB_,opt) // all points
endif
if(reg_type==1)
	return EvSubset_LGR(VB_,opt) // all points
endif
end

function in_train_set(k)
Variable k
Wave /I subset_train
Variable i,n=numpnts(subset_train)
for(i=0;i<n;i+=1)
	if(k==subset_train[i])
		return 1
	endif
endfor
return 0
end

function EvSubset_MLR(v,opt)
Wave v
Variable opt
Variable /G do_regression
Wave /I subset_train,subset_valid,ysubset
Wave yval_set,dval_set
Variable m,n,z=0
Variable ny=numpnts(ysubset)
//
sort v,v
m=purge_repeats(v)
//
	switch(opt)
		case 0: // training set
			if(AutoFill(v,xv,subset_train,0))
				AutoShiftY(yv,ysubset,subset_train,0)
				if(do_regression)
					return LinReg_(xv,yv,yp)
				else
					Wave beta
					return Lin_(xv,yv,beta)
				endif 
			endif
		break
		case 1: // validation set
			if(AutoFill(v,xvb,subset_valid,1))
				AutoShiftY(yvb,ysubset,subset_valid,1)
				return LinRegB_(xvb,yvb,ypb)
			endif
		break
		case -1: // everything
		default:
			n=DimSize(yval_set,1)
			Make /N=(n) /O /I subset_train
			subset_train=p
			if(AutoFill(v,xv,subset_train,0))
				AutoShiftY(yv,ysubset,subset_train,0)
				if(do_regression)
					return LinReg_(xv,yv,yp)
				else
					Wave beta
					return Lin_(xv,yv,beta)
				endif 
			endif
	endswitch
return Inf
end

function AutoShiftY(yv,ysubset,subset,opt)
Wave yv,ysubset,subset
Variable opt
Wave yval_set,yav,ystd
Variable z,z1,i,j,u
Variable ny=numpnts(ysubset)
Variable n=numpnts(subset)
Variable /G flag_shift_y=1
Redimension /N=(ny,n) yv
yv=yval_set[ysubset[p]][subset[q]]
if(!opt)
	Redimension /N=(ny) yav,ystd
	for(j=0;j<ny;j+=1)
	z=0; z1=0
		for(i=0;i<n;i+=1)
			u=yv[j][i]
			z+=u
			z1+=u*u
		endfor
	z/=n
	z1=sqrt(abs(z1/n-z*z))
	yav[j]=z
	ystd[j]=z1
	endfor
endif
yv=(yv-yav[p])/(ystd[p])
end


function renew_best()
Wave VB
Variable /G r2w,rmsw,IT
make_V_comment(VB)
EvVect_(VB,1)
SetYShowProc_(0)
r2vp_to_r2y()
doUpdate
printf ">> iteration %d :: r2=%.3f, rms=%.3f\r\n",IT,1-r2w,rmsw
end



Macro ShowBestSolution()
Make /N=0 /O /R yp0,yv0,ep0
Variable /G reg_type
if(reg_type==0)
	print "Multivariate Multiple Linear Regression"
	hat_eigen_MLR()
endif
if(reg_type==1)
	print "Univariate Multiple Logistic Regression"
	hat_eigen_LGR()
endif
Variable /G nfixed,auto_quad
Variable /G solution_found=1
Variable m=purge_repeats(VB)
if(nfixed)
	printf "%d fixed descriptors\r\n",nfixed
endif
printf "%d random descriptors\r\n",m
EvVect_(VB,0)
EvVect_(VB,1)
Note /K beta
Note beta num2str(reg_type)
make_V_comment(VB)
make_VB_comment()
ShowVB_()
SetYShowProc_(1)
Variable /G rankw=matrixrank(xx,1e10)-1
printf "estimated (weighed) covariance matrix rank is %d\r\n",rankw
Variable /G yshow,pshow=0,cshow=0
Make /N=0 /R /O yc,ycp,xc,yc_,ycp_,xc_
if(!auto_quad)
	//mk_comp_param()
	color_comp_set()
endif
Variable /G rmsw,r2w,rmswb,r2wb
if(reg_type==0)
	printf "MMLR :: R^2(Wilks) = %.4f,  RMS (training)=%.4f, (validation)=%.4f\r\n",1-r2w,rmsw,rmswb
endif
if(reg_type==1)
	printf "LGR :: pseudo-R^2(McFadden) = %.4f, mean deviance (training)=%.4f, (validation)=%.4f\r\n",1-r2w,rmsw,rmswb
endif
//
duplicate /O V_comment_show dexplain_show
dexplain_show=des_explain(V_comment_show)
DoWindow /K Table_coeff
Table_coeff()
if(exists("class_set")==1)
	do_class_set()
endif
end

function fill_yp1(k)
Variable k
Wave /B status_set
Wave yp,ypb,yp1
Wave ysubset
Variable u=WhereInWave(ysubset,k)
if(u<0)
	yp1=Nan
	return 0
endif
Variable i,n=numpnts(status_set)
//
Variable kt=0, kv=0
for(i=0;i<n;i+=1)
	if((status_set[i]==1) && (Check_Nan_y(i)))
		yp1[i]=yp[u][kt]
		kt+=1	
	else
		yp1[i]=ypb[u][kv]
		kv+=1
	endif
endfor
return 1
end

function unshift_y()
Wave yp,ypb,yav,ystd,yv,yvb
Variable /G flag_shift_y
if(flag_shift_y)
	yp=yp*ystd[p]+yav[p]
	ypb=ypb*ystd[p]+yav[p]
	yv=yv*ystd[p]+yav[p]
	yvb=yvb*ystd[p]+yav[p]
	flag_shift_y=0
else
	//print ">> y set is already unshiftted"
endif
end

function SetYShowProc_(opt)
Variable opt
Variable /G beta_sort=opt
	SetYShowProc("",0,"","")
end

function SetYShowProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
PauseUpdate
Silent 1
Variable /G beta_sort
Wave line45,line45_
Wave yval_set
Wave yv,yp,ypb,xv,beta,ystd,yav
Wave /T yname_set
Variable /G yshow
Wave /I ysubset
Variable /G reg_type
Variable m=DimSize(yval_set,1)
Variable n=DimSize(yval_set,0)
Wave /T V_comment_show
Variable /G solution_found
Make /N=(m) /O /R ep1,yp1,yv1,ys1,tv1
ys1=in_train_set(p)
String /G yname_show
Variable i=CheckCircularBelongIndex(ysubset,yshow)
yshow=ysubset[i]
yname_show=yname_set[yshow]
Variable d=DimSize(beta,1)
//
if(reg_type==0)
	Make /O /N=(d)/R beta1
	beta1=beta[i][p]
	unshift_y()
	yv1=yval_set[yshow][p]
	back_y_Calculator(yshow)
	ep1=yv1-yp1
	tv1=label_text(yv1,yp1)
endif
//
if(reg_type==1)
	yshow=0
	yname_show=yname_set[0]
	d-=1
	Make /O /N=(d) /R beta1
	beta1=beta[i][p+1]	
	yv1=yval_set[0][p]
	fill_yp1(0)
	Make /N=(m) /O /R pv1
	pv1=Logit(yp1)
	ep1=Logistic_dev(yv1,yp1)
	tv1=label_text(yv1,pv1)
endif
duplicate /O beta1 beta2
beta2=abs(beta1)
WaveStats /Q yv1
mkline45()
//
if(beta_sort)
	Sort /R beta2,beta2,beta1,V_comment_show
endif
Make /N=(d) /O /T dexplain_show
dexplain_show=des_explain(V_comment_show[p])
Wave r2y
Variable /G r2show=r2y[yshow]
End

function label_text(x,y)
Variable x,y
if(numtype(x))
	return y
endif
return max(x,y)
end

Proc PlotAllMLR_Y()
Variable /G yshow
String /G yname_show
Variable i=0
Variable n=DimSize(ysubset,0)
Silent 1
String s,q
do
	PauseUpdate
	yshow=ysubset[i]
	s=yname_set[yshow]
	yname_show=s
	//print yshow,yname_show
	SetYShowProc("",0,"","")
	s=StringFromList(1,s,"*")
	s=StringFromList(0,s,".")
	sprintf q,"\\Z16\\f01%s (%d)\r\nr^2=%.3f",s,yshow+1,r2y[i]
	TextBox/C/N=text0/F=0/A=MC q
	doUpdate
	SavePICT/O /E=-6/B=288 /P=home as s+".jpg"
i+=1
while(i<n)
end

Proc MovieAllMLR_Y()
Variable /G yshow
String /G yname_show
Variable i=0
Variable n=DimSize(ysubset,0)
Silent 1
String s,q
NewMovie /O/Z /P=home as "MLR.mp4"
do
	PauseUpdate
	yshow=ysubset[i]
	s=yname_set[yshow]
	yname_show=s
	//print yshow,yname_show
	SetYShowProc("",0,"","")
	s=StringFromList(1,s,"*")
	s=StringFromList(0,s,".")
	sprintf q,"\\Z16\\f01%s (%d)\r\nr^2=%.3f",s,yshow+1,r2y[i]
	TextBox/C/N=text0/F=0/A=MC q
	doUpdate
	AddMovieFrame
	i+=1
while(i<n)
CloseMovie
end



function mkline45()
Wave yv1,yp1
Variable z=Inf,z1=-Inf,u,u1
Variable i=0,n=numpnts(yv1)
for(i=0;i<n;i+=1)
	if(!numtype(yv1[i]))
		z=min(z,yv1[i])
		z1=max(z1,yv1[i])	
	endif
	if(!numtype(yp1[i]))
		z=min(z,yp1[i])
		z1=max(z1,yp1[i])	
	endif
endfor
u=(z1+z)/2
u1=abs(z1-z)/2
Make /O line45={u-1.1*u1,u+1.1*u1}
end


Proc ShowVB_()
duplicate /O VB_comment V_comment_show
end


Macro BestVectorFromList()
Variable /G Ng=numpnts(V_comment)
Make /N=(Ng) /O /R VB
VB=WhereInSWave(dname_set,V_comment)
ShowBestSolution()
end

Proc CopyOriginalSet()
SaveSet("copy")
end 

Proc RestoreOriginalSet()
RestoreSet("copy")
String /G y_transform=note(yname_set)
endMacro

Proc RestoreLinSet()
// from unshifted copy
RestoreSet("lin")
endMacro

Proc RestoreQuadSet()
// from unshifted copy
RestoreSet("quad")
endMacro

Proc RestoreFinalSet()
// from unshifted copy
RestoreSet("final")
endMacro


Proc RestoreSet(s)
String s
// from unshifted copy
String list="dname_set;dval_set;yval_set;yname_set;comp_set;smi_set;status_set;av_set;std_set;yav_set;ystd_set"
CopyWaveList(list,s,-1)
Variable /G nfixed=0
mkucomp()
mkparam()
endMacro



proc AcceptLinear()
// saves unshifted copy 
AcceptSet()
SaveSet("lin")
duplicate /O beta beta_set_lin
SaveSolution("lin")
end

Proc AcceptQuad()
// saves unshifted copy
AcceptSet()
SaveSet("quad")
duplicate /O beta beta_set_quad
SaveSolution("quad")
end

Proc AcceptFinal()
// saves unshifted copy
AcceptSet()
SaveSet("final")
duplicate /O beta beta_set_final
SaveSolution("final")
end

function count_independent_des(opt)
Variable opt
Wave /T dname_set
String /G list_des="",list_des_ind=""
String s1,s2,s
Variable k=0
Variable i,m,n=numpnts(dname_set)
for(i=0;i<n;i+=1)
	s=dname_set[i]
	if(strsearch(s,"{",0)>=0)
		k+=1
	endif
	list_des+=des_name_clean(s)+";"
endfor
m=ItemsInList(list_des)
for(i=0;i<m;i+=1)
	s=StringFromList(i,list_des)
	if(strlen(s))
		if(WhichListItem(s,list_des_ind)<0)
			list_des_ind+=s+";"
			if(opt)
				print s
			endif
		endif
	endif
endfor
m=ItemsInList(list_des_ind)
printf ">>> %d descriptors :: %d linear, %d quadratic terms\r\n",m,n-k,k
return m
end

function /S des_name_clean(s)
String s
Variable flag=0
Variable m=strsearch(s,"}^",0)
if(m>0)
	s=s[0,m]
	flag=1
endif
s=replacesymbols(s,"*","")	
s=replacesymbols(s,"}{",";")
s=replacesymbols(s,"{","")
s=replacesymbols(s,"}","")
if(flag)
	s+=";"+s
endif
return s
end


Proc AcceptSet()
String s
PauseUpdate
ShowBestSolution()
Variable m=accept_best()
duplicate /O dval_set_ dval_set
duplicate /O dname_set_ dname_set
duplicate /O xav av_set
duplicate /O xstd std_set
if(reg_type==0)
	duplicate /O yav yav_set
	duplicate /O ystd ystd_set
endif
KillWaves dval_set_,dname_set_
duplicate /O dname_set V_comment
Make /N=(m) /O /R VB
VB=p
Variable /G nfixed=0,Ng=m
Variable /G auto_quad=0
ShowBestSolution()
count_independent_des(0)
end


Proc QuadExpansion()
DoAlert 1,"Accept the set as linear set before expanding it"
DoWindow /K Table_coeff
QuadExpansion_()
simple_stat_test(1)
end

function QuadExpansion_()
Wave dval_set_lin
Wave /T dname_set_lin
Variable m=DimSize(dval_set_lin,0)
Variable n=DimSize(dval_set_lin,1)
Variable m2=m+m*(m+1)/2
Make /N=(m2,n) /O /R dval_set
Make /N=(m2) /O /T dname_set
Variable k,i,j,j1
String s
for(i=0;i<n;i+=1)
	k=m
	for(j=0;j<m;j+=1)
		dval_set[j][i]=dval_set_lin[j][i]
		if(!i)
			s=dname_set_lin[j]
			if(strsearch(s,"*",0)<0)
				s="*"+s
			endif
			dname_set[j]=s
		endif
		for(j1=0;j1<=j;j1+=1)
			dval_set[k][i]=dval_set_lin[j1][i]*dval_set_lin[j][i]
			if(!i)
				if((j==j1))
					sprintf s,"{%s}^2",dname_set_lin[j]
				else
					sprintf s,"{%s}{%s}",dname_set_lin[j1],dname_set_lin[j]
				endif
				dname_set[k]=replacesymbols(s,"*","")
			endif
		k+=1
		endfor
	endfor
endfor
end

Window Table_coeff() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(3,46,675,284) V_comment_show,beta1,dexplain_show
	ModifyTable format(Point)=1,width(Point)=0,alignment(V_comment_show)=0,width(V_comment_show)=117
	ModifyTable rgb(V_comment_show)=(0,0,65280),alignment(beta1)=0,format(beta1)=3,width(beta1)=53
	ModifyTable alignment(dexplain_show)=0,width(dexplain_show)=435,rgb(dexplain_show)=(65280,0,0)
EndMacro


function binary_weight_fixed(p)
Variable p
Wave /B status_set_old
Variable i,n=numpnts(status_set_old)
Variable np=0
Make /N=(n) /R /O color_set
//
duplicate /O status_set_old status_set
status_set=0
for(i=0;i<n;i+=1)
	if(!Check_nan_y(i))
		status_set_old[i]=2
		status_set[i]=2
	else
		if(status_set_old[i]==0)
			status_set[i]=2
		else
			np+=1
		endif
	endif
endfor
np=floor(np*p)
Make /N=(np) /O /I subset_train
random_mask(status_set,subset_train,np,0)
//
Variable m=0
Make /N=0 /O /I subset_valid
for(i=0;i<n;i+=1)
	if((status_set[i]==0) || (status_set_old[i]!=1))
		status_set[i]=0
		Make /N=(m+1) /O /I subset_valid
		subset_valid[m]=i
		m+=1
	endif
endfor
color_set=status_set
end



Window show_MMLR() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(-14,48,581,700) line45 vs line45
	AppendToGraph yv1 vs yp1
	AppendToGraph/L=new ep1 vs yp1
	ModifyGraph mode(yv1)=3,mode(ep1)=1
	ModifyGraph marker(yv1)=8
	ModifyGraph lSize(yv1)=3,lSize(ep1)=2
	ModifyGraph rgb(line45)=(0,0,0),rgb(ep1)=(0,39168,0)
	ModifyGraph msize(yv1)=2
	ModifyGraph opaque(yv1)=1
	ModifyGraph zColor(yv1)={status_set,1,0,RedWhiteBlue}
	ModifyGraph grid(left)=2,grid(bottom)=2
	ModifyGraph zero(new)=4
	ModifyGraph mirror=2
	ModifyGraph minor(left)=1,minor(bottom)=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=30
	ModifyGraph axOffset(bottom)=1.54167
	ModifyGraph axThick=2
	ModifyGraph zeroThick(new)=1
	ModifyGraph lblPos(left)=60,lblPos(new)=63
	ModifyGraph freePos(new)=1
	ModifyGraph axisEnab(left)={0,0.7}
	ModifyGraph axisEnab(new)={0.75,1}
	Label left "exp."
	Label bottom "estimate"
	Label new "res."
	SetAxis/A/N=1/E=2 new
	Cursor/P A yv1 4
	ShowInfo
	ShowTools
	SetVariable setvar0,pos={355.00,627.00},size={90.00,19.00},proc=SetYShowProc,title="show y"
	SetVariable setvar0,value= yshow
	SetVariable setvar1,pos={451.00,627.00},size={90.00,19.00},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={228.00,627.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= yname_show
	Button button0,pos={48.00,625.00},size={60.00,20.00},proc=ButtonProc,title="hist errors"
	Button button0,fSize=11
EndMacro




// ============== export solution as csv file ====================


Proc SaveSolution(suffix)
String suffix
	export_y(suffix)
	export_beta(suffix)
end

function export_beta(c) // needs to be the complete y set
String c
String s
Variable i,j,b,k
Variable /G refnum,refnum1
Variable /G reg_type
Variable opt=reg_type
Wave beta,yav,ystd,ysubset,av_set,std_set
Wave /T yname_set,dname_set
Variable /G rmsw,r2w,rmswb,r2wb
Open /T="TEXT" /P=home refnum as "coeff_"+c+".csv"
Variable n=numpnts(dname_set)
Variable m=numpnts(ysubset)
fprintf refnum,"\r\nCOEFFICIENTS\r\n"
if(reg_type==0)
	fprintf refnum,"MMLR,,,R^2(Wilks),,,RMS(training),,,RMS(validation)"
endif
if(reg_type==1)
	fprintf refnum, "LGR,,,pseudo-R^2(McFadden),,,mean deviance (training),,,mean deviance (validation)"
endif
fprintf refnum,"\r\n,,,%.4f,,,%.4f,,,%.4f",1-r2w,rmsw,rmswb
fprintf refnum,"\r\n\r\ndescriptor"
for(i=0;i<m;i+=1)
	fprintf refnum,",%s",yname_set[ysubset[i]]
endfor
	fprintf refnum,",avg,std,type"
for(j=0;j<n;j+=1) // descriptors
	fprintf refnum,"\r\n%s",dname_set[j]
	for(i=0;i<m;i+=1) // variables		
	      fprintf refnum,",%.6f",beta[i][j+opt]	
	endfor
	s=des_explain(dname_set[j])
	fprintf refnum,",%.6f,%.6f,%s",av_set[j],std_set[j],s
endfor
//
if(opt)
	fprintf refnum,"\r\noffset"
	for(i=0;i<m;i+=1) // variables		
		fprintf refnum,",%.6f",beta[m][0]
	endfor
else
	fprintf refnum,"\r\n=y[av]="
	for(i=0;i<m;i+=1) // variables		
		fprintf refnum,",%.6f",yav[i]
	endfor
	fprintf refnum,"\r\n=y[std]="
	for(i=0;i<m;i+=1) // variables		
		fprintf refnum,",%.6f",ystd[i]
	endfor
endif	
//
Variable u
if(reg_type==0)
fprintf refnum,"\r\n\r\n\r\ndescriptor"
Make /N=(m) /O beta0
beta0=yav
for(i=0;i<m;i+=1)
	fprintf refnum,",%s",yname_set[ysubset[i]]
endfor
	fprintf refnum,",,,type"
for(j=0;j<n;j+=1) // descriptors
	fprintf refnum,"\r\n%s",dname_set[j]
	for(i=0;i<m;i+=1) // variables		
			u=ystd[i]*beta[i][j+opt]/std_set[j]
	      fprintf refnum,",%.6f",u	
	      beta0-=av_set[j]*u
	endfor
	s=des_explain(dname_set[j])
	fprintf refnum,",,,%s",s	
endfor
//
fprintf refnum,"\r\noffset"
if(opt)	
		beta0+=beta[m][0]*ystd
endif	
for(i=0;i<m;i+=1) // variables	
		fprintf refnum,",%.6f",beta0[i]
endfor
endif
Close refnum
end

function export_y(c)  // needs to be the complete y set
String c
Variable /G refnum,yshow
Wave yval_set,yp1,yv
Wave /T yname_set,comp_set,dname_set,smi_set
Wave /B status_set
Wave dval_set,yav,ystd,av_set,std_set
Open /T="TEXT" /P=home refnum as "fit_"+c+".csv"
Variable /G reg_type
String /G nindex
Variable n=numpnts(comp_set)
Variable my=numpnts(yname_set)
Variable m=numpnts(dname_set)
fprintf refnum,"compound,smiles,status"
String s
Variable i,j,b,res,z,zp
for(j=0;j<m;j+=1)
	s=dname_set[j]
	if(strsearch(s,"*",0)>=0)
		fprintf refnum,",%s",s
	endif
endfor
for(j=0;j<my;j+=1)
	if(reg_type==0)
		fprintf refnum,",%s,fit,res",yname_set[j]
	else
		fprintf refnum,",%s,prob,dev,,linear",yname_set[j]
	endif
endfor
//
for(i=0;i<n;i+=1) // compounds
	fprintf refnum,"\r\n%s,%s,%d",comp_set[i],smi_set[i],status_set[i]
	for(j=0;j<m;j+=1) // parameters
		s=dname_set[j]
		if(strsearch(s,"*",0)>=0)
			z=dval_set[j][i]
			fprintf refnum,",%.3f",z
		endif
	endfor
	//
	s=""
	for(j=0;j<my;j+=1) // variables
		if(reg_type==0)
			back_y_Calculator(j)
			z=yval_set[j][i]
			zp=yp1[i]
			sprintf s,"%s,%.3f,%.3f,%.4f",s,z,zp,z-zp
		else
			fill_yp1(j)
			z=Logistic_dev(yval_set[0][i],yp1[i])
			zp=Logit(yp1[i])
			sprintf s,"%s,%.3f,%.3f,%.3f,,%.3f",s,yv[j][i],zp,z,yp1[i]
		endif		
	endfor
	fprintf refnum,"%s",replacesymbols(s,"Nan","")
endfor

// restore currentyp 
if(reg_type==0)
	back_y_Calculator(yshow)
else
	fill_yp1(yshow)
endif	

// selected descriptors
fprintf refnum,"\r\n\r\n\r\ncompound,smiles,status"
for(j=0;j<m;j+=1) // descriptors
	fprintf refnum,",%s",dname_set[j]
endfor
for(i=0;i<n;i+=1) // compounds
	fprintf refnum,"\r\n%s,%s,%d",comp_set[i],smi_set[i],status_set[i]
	for(j=0;j<m;j+=1) // descriptors
		fprintf refnum,",%.6f",dval_set[j][i]
	endfor
endfor
Close refnum
end


n=numpnts(dname_set)
m=numpnts(ysubset)


for(i=0;i<m;i+=1)
	fprintf refnum,",%s",yname_set[ysubset[i]]
endfor
	fprintf refnum,",avg,std,type"

	for(i=0;i<m;i+=1) // variables		
	      fprintf refnum,",%.6f",beta[i][j+opt]	
	endfor
	s=des_explain(dname_set[j])
	fprintf refnum,",%.6f,%.6f,%s",av_set[j],std_set[j],s
endfor


//==============================


// expects <name><separator><smi> format

Proc ExportCompList(sd,sy)
String sd=""
String sy=""
Prompt sd,"descriptors"
Prompt sy,"y's"
Variable /G refnum
Open /P=home refnum as "smiles_list.csv"
Variable i,n
export_smiles_des(sd,sy)
Close /A
end

function export_smiles_des(sd,sy)
String sd,sy
Variable /G refnum
Wave /T comp_set,smi_set,dname_set,yname_set
Wave dval_set,yval_set,status_set
Variable u,i,n,nd=0,ny=0
String s

fprintf refnum,"name,smiles,status"

n=ItemsInList(sd)
for(i=0;i<n;i+=1)
	u=find_des_index(StringFromList(i,sd))
	if(u>=0)
		Make /N=(nd+1) /O /I indd
		indd[nd]=u
		nd+=1
		fprintf refnum,",%s",dname_set[u]
	endif
endfor

n=ItemsInList(sy)
for(i=0;i<n;i+=1)
	u=find_y_index(StringFromList(i,sy))
	if(u>=0)
		Make /N=(ny+1) /O /I indy
		indy[ny]=u
		ny+=1
		fprintf refnum,",%s",yname_set[u]
	endif
endfor

n=numpnts(comp_set)
for(i=0;i<n;i+=1)
	s=smi_set[i]
	if(strlen(s))
		fprintf refnum,"\r\n%s,%s,%d",comp_set[i],s,status_set[i]
		for(u=0;u<nd;u+=1)
			fprintf refnum,",%g",dval_set[indd[u]][i]
		endfor
		for(u=0;u<ny;u+=1)
			fprintf refnum,",%g",yval_set[indy[u]][i]
		endfor
	endif
endfor

KillWaves /Z indd,indy
end


Proc MkUniqueCompList()
Silent 1
PauseUpdate
LoadWave/J/N=wave/O/K=2/Q/V={","," $",0,1}
Variable i=1,n=numpnts(wave0),m=1
Make /N=1 /T /O comp_set
comp_set[0]=wave0[0]
if(exists("wave1")==1)
	Make /N=1 /T /O smi_set
	smi_set[0]=wave1[0]
endif
String s
do
	s=replacesymbols(wave0[i]," ","")
	if(WhereInSWave(comp_set,s)<0)
		Make /N=(m+1) /T /O comp_set,smi_set
		comp_set[m]=s
		if(exists("wave1")==1)
			smi_set[m]=wave1[i]
		endif
		m+=1
	endif
i+=1
while(i<n)
KillWaves wave0
if(exists("wave1")==1)
	KillWaves wave1
	Duplicate /O comp_set,smi_set
	smi_set=""
endif
printf ">>found %d unique compounds\r\n",m
if(n!=m)
	Open /T="TEXT" refnum as "parsed_comp_list.csv"
	i=0
	do
		fprintf refnum,"%s,%s\r\n",comp_set[i],smi_set[i]
	i+=1
	while(i<m)
	Close refnum
endif
end

function check_std_set()
Wave /T dname_set
Wave std_set
Variable i
Variable n=numpnts(std_set)
for(i=0;i<n;i+=1)
	if(std_set[i]<1e-6)
		dname_set[i]="#"+dname_set[i]
	endif
endfor
end


Macro ColorizeComp()
color_comp_set()
DoWindow /K show_qspr_shifted
DoWindow /K show_qspr_shifted_color
show_qspr_shifted_color()
end

function color_comp_set()
String  /G ucomp_set
Wave /T comp_set
Variable n=numpnts(comp_set)
Make /N=(n) /O /R color_set
//color_set=WhichListItem(comp_set[p],ucomp_set)
end

Window show_MMLR_color() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(8,48,578,584) line45 vs line45
	AppendToGraph yv1 vs yp1
	ModifyGraph mode(yv1)=3
	ModifyGraph marker(yv1)=19
	ModifyGraph lSize(yv1)=3
	ModifyGraph rgb(line45)=(0,0,0)
	ModifyGraph msize(yv1)=5
	ModifyGraph opaque(yv1)=1
	ModifyGraph zColor(yv1)={color_set,*,*,Rainbow}
	ModifyGraph grid=2
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=25
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph lblLatPos(bottom)=-1
	Label left "exp. (shifted)"
	Label bottom "estimate (shifted)"
	Cursor/P A yv1 1
	ShowInfo
	SetVariable setvar0,pos={364.00,514.00},size={90.00,19.00},proc=SetYShowProc,title="show y"
	SetVariable setvar0,value= yshow
	SetVariable setvar1,pos={460.00,514.00},size={90.00,19.00},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={237.00,514.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= yname_show
EndMacro

Window bar_plot_MMLR() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(10,46,585,670) yv1,yp1
	AppendToGraph/L=new ep1
	ModifyGraph mode(yv1)=1,mode(yp1)=3,mode(ep1)=1
	ModifyGraph marker(yv1)=8,marker(yp1)=19
	ModifyGraph lSize(yv1)=3,lSize(ep1)=2
	ModifyGraph rgb(ep1)=(0,39168,0)
	ModifyGraph msize(yv1)=2
	ModifyGraph opaque(yv1)=1
	ModifyGraph zColor(yp1)={ys1,-1,1,BlueHot}
	ModifyGraph grid(left)=2,grid(bottom)=2
	ModifyGraph zero(new)=4
	ModifyGraph mirror(left)=2,mirror(bottom)=2
	ModifyGraph minor(left)=1,minor(bottom)=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=28
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph zeroThick(new)=1
	ModifyGraph lblPos(left)=69,lblPos(new)=65
	ModifyGraph lblLatPos(bottom)=-10,lblLatPos(new)=-5
	ModifyGraph freePos(new)=2
	ModifyGraph axisEnab(left)={0,0.7}
	ModifyGraph axisEnab(new)={0.75,1}
	Label left "y's"
	Label bottom "compounds"
	Label new "res."
	SetAxis/A/N=1 left
	SetAxis/A/N=1 bottom
	SetAxis/A/N=1/E=2 new
	Cursor/P A yp1 7
	ShowInfo
	Legend/C/N=text0/J/A=MC/X=18.63/Y=27.88 "\\Z09\rblue - training\ryellow - validation\r\\s(yv1) exp"
	SetVariable setvar1,pos={336.00,602.00},size={90.00,19.00},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={207.00,602.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= yname_show
	Button button1,pos={499.00,600.00},size={25.00,20.00},proc=ShiftRight,title=">>"
	Button button2,pos={461.00,600.00},size={25.00,20.00},proc=ShiftLeft,title="<<"
	SetVariable setvar0,pos={110.00,603.00},size={90.00,19.00},proc=SetYShowProc,title="show y"
	SetVariable setvar0,value= yshow
EndMacro

Window bar_plot_MMLR_names() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(10,46,1180,670) yv1,yp1,tv1
	ModifyGraph mode(yv1)=1,mode(yp1)=3,mode(tv1)=3
	ModifyGraph marker(yv1)=8,marker(yp1)=19
	ModifyGraph lSize(yv1)=3
	ModifyGraph rgb(tv1)=(0,0,0)
	ModifyGraph msize(yv1)=2,msize(tv1)=4
	ModifyGraph opaque(yv1)=1
	ModifyGraph zColor(yp1)={ys1,-1,1,BlueHot}
	ModifyGraph textMarker(tv1)={comp_set,"default",0,90,1,0.00,50.00}
	ModifyGraph grid=2
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=28
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph lblPos(left)=69
	ModifyGraph lblLatPos(bottom)=-10
	ModifyGraph axisEnab(left)={0,0.9}
	Label left "y's"
	Label bottom "compounds"
	SetAxis/A/N=1 left
	SetAxis/N=1 bottom -1,64
	Cursor/P A yp1 7
	ShowInfo
	TextBox/C/N=text10/A=MC/X=-27.65/Y=45.97 "\\Z09\rblue - training\ryellow - validation\r\\s(yv1) exp"
	SetVariable setvar1,pos={336.00,602.00},size={90.00,19.00},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={207.00,602.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= yname_show
	Button button1,pos={499.00,600.00},size={25.00,20.00},proc=ShiftRight,title=">>"
	Button button2,pos={461.00,600.00},size={25.00,20.00},proc=ShiftLeft,title="<<"
	SetVariable setvar0,pos={110.00,603.00},size={90.00,19.00},proc=SetYShowProc,title="show y"
	SetVariable setvar0,value= yshow
EndMacro

Window bar_plot_LGR() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(10,46,615,671) yv1,pv1
	AppendToGraph/L=new ep1
	ModifyGraph mode(yv1)=1,mode(pv1)=3,mode(ep1)=1
	ModifyGraph marker(yv1)=8,marker(pv1)=19
	ModifyGraph lSize(yv1)=3,lSize(ep1)=2
	ModifyGraph rgb(pv1)=(0,0,65280),rgb(ep1)=(0,39168,0)
	ModifyGraph msize(yv1)=2,msize(pv1)=3
	ModifyGraph opaque(yv1)=1
	ModifyGraph zColor(pv1)={ys1,-1,1,BlueHot}
	ModifyGraph grid(left)=2,grid(bottom)=2
	ModifyGraph zero(new)=4
	ModifyGraph mirror(left)=2,mirror(bottom)=2
	ModifyGraph minor(left)=1,minor(bottom)=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=28
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph zeroThick(new)=1
	ModifyGraph lblPos(left)=69,lblPos(new)=65
	ModifyGraph lblLatPos(bottom)=-10,lblLatPos(new)=-5
	ModifyGraph freePos(new)=2
	ModifyGraph axisEnab(left)={0,0.7}
	ModifyGraph axisEnab(new)={0.75,1}
	Label left "y's and p's"
	Label bottom "compounds"
	Label new "deviance"
	SetAxis left 0,1.05
	SetAxis bottom -1,50
	SetAxis/A/N=1/E=2 new
	Cursor/P A pv1 17
	ShowInfo
	Legend/C/N=text0/J/A=MC/X=47.42/Y=36.59 "\\Z09\\s(pv1) probability \rblue - training\ryellow - validation\r\\s(yv1) binary y\r\\s(ep1) deviance"
	SetVariable setvar1,pos={336.00,602.00},size={90.00,19.00},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={207.00,602.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= yname_show
	Button button1,pos={499.00,600.00},size={25.00,20.00},proc=ShiftRight,title=">>"
	Button button2,pos={461.00,600.00},size={25.00,20.00},proc=ShiftLeft,title="<<"
EndMacro


Window bar_plot_LGR_names() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(10,46,615,671) yv1,pv1,tv1
	ModifyGraph mode(yv1)=1,mode(pv1)=3,mode(tv1)=3
	ModifyGraph marker(yv1)=8,marker(pv1)=19
	ModifyGraph lSize(yv1)=3
	ModifyGraph rgb(pv1)=(0,0,65280),rgb(tv1)=(0,0,0)
	ModifyGraph msize(yv1)=2,msize(pv1)=3,msize(tv1)=4
	ModifyGraph opaque(yv1)=1
	ModifyGraph zColor(pv1)={ys1,-1,1,BlueHot}
	ModifyGraph textMarker(tv1)={comp_set,"default",0,90,1,0.00,10.00}
	ModifyGraph grid(left)=2
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=28
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph lblPos(left)=69
	ModifyGraph lblLatPos(bottom)=-10
	ModifyGraph axisEnab(left)={0,0.85}
	Label left "y's and p's"
	Label bottom "compounds"
	SetAxis left 0,1
	SetAxis bottom -1,20
	ShowInfo
	Legend/C/N=text0/J/A=MC/X=18.24/Y=49.22 "\\Z09\\s(pv1) probability \rblue - training\ryellow - validation\r\\s(yv1) binary y"
	SetVariable setvar1,pos={336.00,602.00},size={90.00,19.00},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={207.00,602.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= yname_show
	Button button1,pos={499.00,600.00},size={25.00,20.00},proc=ShiftRight,title=">>"
	Button button2,pos={461.00,600.00},size={25.00,20.00},proc=ShiftLeft,title="<<"
EndMacro


Proc ShiftRight(ctrlName) : ButtonControl
	String ctrlName
ShiftAxis(1)
End

Proc ShiftLeft(ctrlName) : ButtonControl
	String ctrlName
ShiftAxis(-1)
End

Proc ShiftAxis(z)
Variable z
Variable /G V_min,V_max
GetAxis /Q bottom
Variable /G xcent=(V_max+V_min)/2
Variable /G xsw=abs(V_max-V_min)/2
	xcent+=xsw*z
	SetAxis bottom xcent-xsw,xcent+xsw
end



Proc HistError(opt)
Variable opt
Variable /G V_rms
PauseUpdate
Duplicate /O yp1 ep1
ep1=yp1-yv1
WaveStats /Q ep1
if(opt)
	ep1/=V_rms
endif
Variable n=numpnts(ep1)
Make /O /R /N=(n/10) eh1
Histogram/B=1 ep1,eh1
eh1*=100/n
DoWindow /K hist_err_shifted
hist_err_shifted()
if(opt)
	Label bottom "error/rms"
else
	Label bottom "error (shifted)"
endif
end


Window hist_err_shifted() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,42.5,312.75,367.25) eh1
	ModifyGraph mode=5
	ModifyGraph zero(bottom)=3
	ModifyGraph mirror=2
	ModifyGraph minor(bottom)=1
	Label bottom "error/rms"
	SetAxis/A/N=1/E=1 left
	SetAxis/A/N=1/E=2 bottom
	ShowInfo
EndMacro

Proc ButtonProc(ctrlName) : ButtonControl
	String ctrlName
HistError(1)
End


Proc TestDim(m)
Variable m=2000
Prompt m,"number of iterations"
Silent 1
PauseUpdate
Make /O /R test_par={5,10,15,20,25,30,40,50,60}
Variable n=numpnts(test_par)
Variable ny=DimSize(yval_set,0)
Variable /G Ng,NT=m
Make /N=(n) /O /R rms_par,r2w_par,rank_par,ng_par
rms_par=nan
r2w_par=nan
rank_par=nan
ng_par=nan
rank_par[0]=0
ng_par[0]=0
Prepare_MRL_genetics(); NT=m
Variable i=0
Variable /G rmsw=nan,r2w=nan,rankw=0
do
	DoWindow /K Table_coeff
	Ng=test_par[i]
	printf "================   %d parameters =================\r\n",Ng	
	make_random_vector()
	Genetic(m)
	ShowBestSolution()
	rms_par[i]=rmsw
	r2w_par[i]=r2w
	rank_par[i]=rankw
	ng_par[i]=ng
i+=1
while(i<n)
end

Proc TestDimLQ(m)
Variable m=2000
Prompt m,"number of iterations"
Silent 1
PauseUpdate
Make /O /R test_par={0,10,15,20,25,30,40,50,60}
Variable n=numpnts(test_par)
Variable ny=DimSize(yval_set,0)
Make /N=(n) /O /R rms_par_lin,r2w_par_lin,rank_par_lin,ng_par_lin
Make /N=(n) /O /R rms_par_quad,r2w_par_quad,rank_par_quad,ng_par_quad
rms_par_lin=nan
r2w_par_lin=nan
rank_par_lin=nan
ng_par_lin=nan
rank_par_lin[0]=0
ng_par_lin[0]=0
//
rms_par_quad=nan
r2w_par_quad=nan
rank_par_quad=nan
ng_par_quad=nan
rank_par_quad[0]=0
ng_par_quad[0]=0
//
Variable /G Ng,NT
Variable i=1,k
Variable /G rmsw,r2w,rankw
do
		Ng=test_par[i]
		printf "================   %d parameters =================\r\n",Ng	
		DoWindow /K Table_coeff
		RestoreOriginalSet()
		Prepare_MRL_genetics(); NT=m
		Ng=test_par[i]
		make_random_vector()
		Genetic(m)
		//
		AcceptLinear()
		rms_par_lin[i]=rmsw
		r2w_par_lin[i]=r2w
		rank_par_lin[i]=rankw
		ng_par_lin[i]=ng
		//
		QuadExpansion()
		UnstarAllDescriptors()
		k=0
		//do
			Prepare_MRL_genetics(); NT=m
			if(!k)
				Ng=2*test_par[i]
			else
				Ng=floor(rankw*.9)
			endif
			Genetic(m)
			ShowBestSolution()
			if(rankw==Ng)
				//break
			endif		
		//k+=1	
		//while(1)
		AcceptQuad()
		rms_par_quad[i]=rmsw
		r2w_par_quad[i]=r2w
		rank_par_quad[i]=rankw
		ng_par_quad[i]=ng
i+=1
while(i<n)
end


Window r2_rms_vs_par() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,42.5,432,469.25) rms_par vs ng_par
	AppendToGraph/R r2w_par vs test_par
	ModifyGraph mode=4
	ModifyGraph marker(rms_par)=8,marker(r2w_par)=5
	ModifyGraph lSize=3
	ModifyGraph rgb(r2w_par)=(0,0,65280)
	ModifyGraph msize=6
	ModifyGraph mrkThick=2
	ModifyGraph opaque=1
	ModifyGraph mirror(bottom)=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axOffset(bottom)=0.583333
	ModifyGraph axThick=2
	ModifyGraph axRGB(left)=(65280,0,0),axRGB(right)=(0,0,65280)
	ModifyGraph tlblRGB(left)=(65280,0,0),tlblRGB(right)=(0,0,65280)
	ModifyGraph alblRGB(left)=(65280,0,0),alblRGB(right)=(0,0,65280)
	Label left "RMSE (shifted)"
	Label bottom "# of descriptors"
	Label right "1-R\\S2\\M(Wilks)"
	SetAxis/A/E=1 left
	SetAxis/A/N=1/E=1 bottom
	SetAxis/A/E=1 right
	Cursor/P A rms_par 3
	ShowInfo
EndMacro


Window rank_vs_par() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,42.5,432,469.25) rank_par vs ng_par
	ModifyGraph mode=4
	ModifyGraph marker=8
	ModifyGraph lSize=3
	ModifyGraph rgb=(0,52224,0)
	ModifyGraph msize=6
	ModifyGraph mrkThick=2
	ModifyGraph opaque=1
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axOffset(bottom)=0.583333
	ModifyGraph axThick=2
	Label left "rank of covariance matrix"
	Label bottom "# of descriptors"
	SetAxis/A/E=1 left
	SetAxis/A/N=1/E=1 bottom
	ShowInfo
EndMacro



Window r2_rms_vs_par_quad() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,42.5,432,469.25) rms_par_lin vs ng_par_lin
	AppendToGraph/R r2w_par_lin vs ng_par_lin
	AppendToGraph rms_par_quad vs ng_par_lin
	AppendToGraph/R r2w_par_quad vs ng_par_lin
	ModifyGraph mode=4
	ModifyGraph marker(rms_par_lin)=8,marker(r2w_par_lin)=5,marker(rms_par_quad)=19
	ModifyGraph marker(r2w_par_quad)=16
	ModifyGraph lSize=3
	ModifyGraph lStyle(rms_par_lin)=3,lStyle(r2w_par_lin)=3
	ModifyGraph rgb(r2w_par_lin)=(0,0,65280),rgb(r2w_par_quad)=(0,0,65280)
	ModifyGraph msize=6
	ModifyGraph mrkThick(rms_par_lin)=2
	ModifyGraph opaque(rms_par_lin)=1,opaque(r2w_par_lin)=1
	ModifyGraph mirror(bottom)=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axOffset(bottom)=0.583333
	ModifyGraph axThick=2
	ModifyGraph axRGB(left)=(65280,0,0),axRGB(right)=(0,0,65280)
	ModifyGraph tlblRGB(left)=(65280,0,0),tlblRGB(right)=(0,0,65280)
	ModifyGraph alblRGB(left)=(65280,0,0),alblRGB(right)=(0,0,65280)
	Label left "RMSE (shifted)"
	Label bottom "# of descriptors"
	Label right "1-R\\S2\\M(Wilks)"
	SetAxis/A/E=1 left
	SetAxis/A/N=1/E=1 bottom
	SetAxis/A/E=1 right
	ShowInfo
	TextBox/N=text0/F=0/A=MC/X=-33.99/Y=-40.09 "\\Z09\\s(rms_par_lin) rms, lin\r\\s(r2w_par_lin) r2, lin"
	AppendText "\\s(rms_par_quad) rms, quad\r\\s(r2w_par_quad) r2, quad"
EndMacro

Window rank_vs_par_quad() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(5.25,42.5,432,469.25) rank_par_lin,rank_par_quad vs ng_par_lin
	ModifyGraph mode=4
	ModifyGraph marker(rank_par_lin)=8,marker(rank_par_quad)=19
	ModifyGraph lSize=3
	ModifyGraph lStyle(rank_par_lin)=2
	ModifyGraph rgb=(0,52224,0)
	ModifyGraph msize=6
	ModifyGraph mrkThick=2
	ModifyGraph opaque=1
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph axOffset(bottom)=0.583333
	ModifyGraph axThick=2
	Label left "rank of covariance matrix"
	Label bottom "# of descriptors"
	SetAxis/A/E=1 left
	SetAxis/A/N=1/E=1 bottom
	ShowInfo
	TextBox/N=text0/F=0/A=MC/X=35.31/Y=-43.63 "\\Z08\\s(rank_par_lin) rank, lin\r\\s(rank_par_quad) rank, quad"
EndMacro



Window Table_dim_par() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,42.5,212.25,218) test_par,rms_par,r2w_par,rank_par
	ModifyTable width(Point)=20,width(test_par)=27,format(rms_par)=3,digits(rms_par)=4
	ModifyTable width(rms_par)=39,format(r2w_par)=3,digits(r2w_par)=4,width(r2w_par)=45
	ModifyTable width(rank_par)=48,rgb(rank_par)=(0,0,52224)
EndMacro

Window Table_dim_par_quad() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,42.5,461.25,218) test_par,rms_par_lin,r2w_par_lin,ng_par_lin,rank_par_lin
	AppendToTable rms_par_quad,r2w_par_quad,ng_par_quad,rank_par_quad
	ModifyTable width(Point)=20,width(test_par)=47,format(rms_par_lin)=3,digits(rms_par_lin)=4
	ModifyTable width(rms_par_lin)=51,format(r2w_par_lin)=3,digits(r2w_par_lin)=4,width(r2w_par_lin)=56
	ModifyTable width(ng_par_lin)=33,rgb(ng_par_lin)=(65280,0,0),width(rank_par_lin)=33
	ModifyTable rgb(rank_par_lin)=(0,0,52224),format(rms_par_quad)=3,digits(rms_par_quad)=4
	ModifyTable width(rms_par_quad)=54,format(r2w_par_quad)=3,digits(r2w_par_quad)=4
	ModifyTable width(r2w_par_quad)=53,width(ng_par_quad)=35,rgb(ng_par_quad)=(65280,0,0)
	ModifyTable width(rank_par_quad)=38,rgb(rank_par_quad)=(0,0,52224)
EndMacro



Window Table_mordered() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5,43,603,745) dname_m,type_m,d_m
	ModifyTable format(Point)=1,width(Point)=30,size(dname_m)=9,alignment(dname_m)=0
	ModifyTable width(dname_m)=93,size(type_m)=9,alignment(type_m)=0,width(type_m)=50
	ModifyTable rgb(type_m)=(65280,0,0),size(d_m)=9,alignment(d_m)=0,width(d_m)=479
EndMacro

function /S des_explain(s)
String s
String q=des_name_clean(s)
if(!strlen(q))
	return ""
endif
//
String s1=StringFromList(0,q)
String s2=StringFromList(1,q)
if(strsearch(s,"{",0)>=0)
	if(strsearch(s,"}^",0)>=0)
		sprintf s,"{%s}^2",des_explain_(s1)
	else
		sprintf s,"{%s}{%s}",des_explain_(s1),des_explain_(s2)
	endif
	return s
else
	return des_explain_(s)
endif
end

function /S des_explain_(s)
String s
Wave /T dname_m,type_m,d_m
Variable m,a,b,c
String q
s=replacesymbols(s," ","")
s=replacesymbols(s,"(+)","")
s=replacesymbols(s,"(-)","")
if(strsearch(s,".G98",0)>=0)
	return "Gaussian"
endif
m=WhereInSWave(dname_m,s)
if(m>=0)
	sprintf s,"%s: %s",type_m[m],d_m[m]
	return s
endif
if(!strsearch(s,"Ring",0))
	sscanf s,"Ring%d%c%d",a,b,c
	if(c)
		sprintf s,"1D: %d-atom rings containing %d %c-atoms",a,c,b
	else
		if(b)
			sprintf s,"1D: %d-atom heterocycles",a
		else
			sprintf s,"1D: %d-atom rings",a
		endif
	endif
	return s
endif
return ""
end





Macro DoAllYPlots()
Variable ny=numpnts(ysubset)
String /G yname_show
Variable /G yshow
Variable i=0
String s,sp
do
	yshow=ysubset[i]
	SetYShowProc_(0)
	DoUpdate
	s=replacesymbols(yname_show,"*","")
	sp="^"+s
	duplicate /O yv1_   $s
	duplicate /O yp1_	  $sp
i+=1
while(i<ny)
end



function eval_from_param(k,ky,wname,wval) // uses final rather than current sets
Variable k,ky
Wave /T wname
Wave wval
Wave /T dname_set_final
Wave av_set_final,std_set_final,beta_final
Wave beta_set_final
Variable /G reg_type
Variable opt=reg_type
Variable j,z,u,z1
String s
Variable m=numpnts(dname_set_final)
u=0
for(j=0;j<m;j+=1)
	s=dname_set_final[j]
	z=param2val(k,s,wname,wval)
	z1=(z-av_set_final[j])/std_set_final[j]
	u+=beta_set_final[ky][j+opt]*z1
	//print j,s,z,z1,av_set_final[j],std_set_final[j],beta_set_final[ky][j+opt],u
endfor
if(reg_type==1)
	u+=beta_set_final[ky][0]
	u=1/(1+exp(-u))
endif	
return u
end

function param2val(k,s,wname,wval)
Variable k //compound
String s // parameter name
Wave /T wname // parameter names
Wave wval // parameter wave
String  s1,s2
Variable z1,z2
s=des_name_clean(s)
s1=StringFromList(0,s)
s2=StringFromList(1,s)
if(!strlen(s2))
	return identify_param(k,s,wname,wval)
else
	return identify_param(k,s1,wname,wval)*identify_param(k,s2,wname,wval)	
endif
end

function identify_param(k,s,wname,wval)
Variable k
String s
Wave /T wname  // parameter names
Wave wval    // parameters
Wave /T pname
Wave pval
Variable i,n=numpnts(pname)
if(n)
	for(i=0;i<n;i+=1)
	if(stringmatch(s,pname[i]))
		return pval[i]
	endif
endfor
endif
//
i=WhereInSWave(wname,s)
if(i>=0)
	return wval[i][k]
else
	return nan
endif
end


//=============================================================

Proc EvalLoadedSet(type)
Variable type=0
Prompt type,"regression type: 0 - MMLR, 1 - logistic"
Variable /G reg_type=type
DoAlert 0,"final set is used as solution"
Variable /G refout
Variable /G n_test=0
Open /P=home refout as "eval_set.csv"
duplicate /O yname_set_final yname_set
	if(eval_from_current_param(1)==0)
		eval_from_current_param(0)
	endif
Close(refout)
end

Proc EvalLoadedDesDir(type)
Variable type=0
Prompt type,"regression type: 0 - MMLR, 1 - logistic"
Variable /G reg_type=type
DoAlert 0,"final set is used as solution"
Silent 1
PauseUpdate
Variable /G n_test=0
Variable /G refnum,refout
String /G s_filename,s_path
String s
Variable m
Open /P=home refout as "eval_set.csv"
Open /R /P=home /T=".des" refnum
Close(refnum)
removefilepath_(s_filename,":")
Newpath /O des s_path
String /G directory_list=IndexedFile(des,-1,".des")
Variable i=0
eval_from_current_param(1)
String /G smiles=""
Variable /G table_inp=0
Make /N=1 /O /R status_set={1}
do
	s=StringFromList(i,directory_list)
	if(strlen(s))
		Open /R /P=des refnum as s
		s= removefileextension(s)
		s=replacesymbols(s,"_out","")
		printf ">>%s\r\n",s
		Load_descriptors()
		Close(refnum)
		m=numpnts(dname)
		dname=replacesymbols(dname," ","")
		duplicate /O dname dname_set
		Make /N=(m,1) /O /R dval_set
		Make /N=1 /O /T comp_set,smi_set
		dval_set[][0]=dval[p]
		comp_set[0]=s
		smi_set[0]=smiles
		mkucomp()
		duplicate /O yname_set_final yname_set
		eval_from_current_param(0)
	else
		break
	endif
i+=1
while(1)
Close(refout)
end

function eval_from_current_param(opt)
Variable opt
Variable /G refout
Wave /T yname_set_final
Variable i,j,k,z,score=0,j1
String s
Variable /G n_test
Execute "UnstarAllDescriptors()"
Execute "UnhashAllDescriptors()"
choose_y(yname_set_final)
Make /N=0 /T /O pname
Make /N=0 /O pval
Wave /T dname_set_final,yname_set_final
Wave yav_set_final,ystd_set_final
Variable n=DimSize(dval_set,1)
Variable m=numpnts(dname_set)
Variable m1=numpnts(dname_set_final)
Variable ny=numpnts(ysubset)
Make /N=(m1) /O /R locate_set
locate_set=nan
Wave /T dname_set,dname_set_final,yname_set_final
Wave /T comp_set,comp_set_final
Wave dval_set,yval_set_final
//
if(opt)
	fprintf refout,"compound,smiles,"
endif
for(j=0;j<m1;j+=1)
	s=dname_set_final[j]
	s=replacesymbols(s,"*","")
	s=replacesymbols(s,"#","")
	if(strsearch(s,"{",0)<0)
		k=WhereInSWave(dname_set,s)
		locate_set[j]=k
		if(k<0)
			printf "==> %s (final) is not found\r\n",s
			score+=1
		endif		
		if(opt)
			fprintf refout,"%s,",s
		endif
	endif
endfor
if(score)
	return score
endif
//
Wave ysubset
if(opt)
for(j=0;j<ny;j+=1)
	s=yname_set_final[ysubset[j]]
	s=removefileextension(s)
	fprintf refout,"%s.y%d,%s.p%d,",s,j+1,s,j+1
endfor
endif
//
Wave /T smi_set
if(!opt) // evaluate
for(i=0;i<n;i+=1)
	fprintf refout,"\r\n%s,%s,",comp_set[i],smi_set[i]
	for(j=0;j<m1;j+=1)
		k=locate_set[j]
		if(k>=0)
			fprintf refout,"%.4f,",dval_set[k][i]
		endif
		if(k==-1)
			fprintf refout,","
		endif
	endfor	
	//
	Make /N=(n_test+1) /T /O comp_test,smi_test
	Make /N=((n_test+1)*ny) /R /O p_test
	comp_test[n_test]=comp_set[i]
	smi_test[n_test]=smi_set[i]
	//
	Variable /G reg_type
	k=WhereInSWave(comp_set_final,comp_set[i])
	for(j=0;j<ny;j+=1)
		j1=ysubset[j]
		z=eval_from_param(i,j,dname_set,dval_set)
		if(reg_type==0)
			z=z*ystd_set_final[j]+yav_set_final[j]
		endif
		if(k>=0)
			sprintf s,"%.3f,%.3f,",yval_set_final[j1][k],z
		else
			sprintf s,",%.3f,",z
		endif
		fprintf refout,replacesymbols(s,"Nan","")
		p_test[n_test*ny+j]=z
	endfor
	//
	n_test+=1
endfor
endif
//
Variable /G yshow
Redimension /N=(n_test,ny) p_test
MatrixTranspose p_test
Make /N=(n_test) /O /R p_test_show
p_test_show=p_test[yshow][p]
return 0
end

Proc yCalculator(i,fun)
Variable i=1,opt=1
String fun="log"
Prompt i,"which y (1,..)"
Prompt fun,"which function"
// edit this string to change dependency patterns
Silent 1
PauseUpdate
String /G y_transform
i-=1
String s=yname_set[i]
s=replacesymbols(s,"*","")
s=removefileextension(s)
sprintf s,"*%s(%s).y%d",fun,s,ny+1

Variable ny=DimSize(yval_set,0)
Variable n=DimSize(yval_set,1)
Make /O /T /N=(ny+1) yname_set
yname_set[ny]=s
Make /O /R /N=(ny+1,n) yval_set_
Variable j=0
do
	yval_set_[j][]=yval_set[j][q]
j+=1
while(j<ny)
yval_set_[ny][]=y_fun(yval_set[i][q],fun,1)
y_dependencies_remove(i)
sprintf s,"%d;%d;%s,",i,ny,fun
y_transform+=s
duplicate /O yval_set_ yval_set
KillWaves yval_set_
Note /K yname_set
Note yname_set, y_transform
end

function y_dependencies_remove(i)  // remove all dependencies of y[i+1]
Variable i
String /G y_transform
Variable k,n=ItemsInList(y_transform,",")
String s,s1=""
for(k=0;k<n;k+=1)
	s=StringFromList(k,y_transform,",")
	k0=GetFromList(0,s)
		if(k0!=i)
			s1+=s+","
		endif
	endfor
y_transform=s1
end

Proc yDelete(i)
Variable i=1
Prompt i,"which y (1,..)"
Silent 1
PauseUpdate
String /G y_transform
contract_y(i-1)
duplicate /O yname_set_,yname_set
duplicate /O yval_set_,yval_set
KillWaves yval_set_,yname_set_
Note /K yname_set
Note yname_set y_transform
end

function contract_y(i)
Variable i
Wave yval_set
Wave /T yname_set
Variable ny=DimSize(yval_set,0)
Variable n=DimSize(yval_set,1)
Make /N=(ny-1) /T yname_set_
Make /N=(ny-1,n)  /O yval_set_
Variable j,m=0
for(j=0;j<ny;j+=1)
	if(j!=i)
		yval_set_[m][]=yval_set[j][q]
		yname_set_[m]=yname_set[j]
		m+=1
	endif
endfor
//
String /G y_transform
String s,s1=""
Variable k,k0,k1
n=ItemsInList(y_transform,",")	
for(k=0;k<n;k+=1)
	s=StringFromList(k,y_transform,",")
	k0=GetFromList(0,s)
	k1=GetFromList(1,s)
	if((k!=i) && (k1!=i))
		s1+=s+","
	endif
endfor
y_transform=s1
end

function back_y_Calculator(ky)  // always call for shifted y
Variable ky
Wave yp1
Wave /T yname_set
String s
Variable k0,k1,k,m
Variable j,n=DimSize(yval_set,1)
String /G y_transform=note(yname_set)
m=ItemsInList(y_transform,",")
Variable /G reg_type
if(reg_type==0)
for(k=0;k<m;k+=1)
	s=StringFromList(k,y_transform,",")
	k0=GetFromList(0,s)
	k1=GetFromList(1,s)
	s=StringFromList(2,s)
	// back
	if(k0==ky)				
		if( (strsearch(yname_set[k0],"*",0)<0) && (strsearch(yname_set[k1],"*",0)>=0) )
		      fill_yp1(k1)
		      yp1=y_fun(yp1,s,-1)
		      return 1
		endif
	endif
	//
	if(k1==ky)	
		if( (strsearch(yname_set[k0],"*",0)>=0) && (strsearch(yname_set[k1],"*",0)<0) )
		      fill_yp1(k0)
		      yp1=y_fun(yp1,s,+1)
		      return 1
		endif
	endif
endfor
endif
fill_yp1(ky)
return 0
end


function y_fun(x,fun,d)
Variable x,d
String fun
Variable z=nan
if(d>0)
	strswitch(fun)
		case "log":
			z=log(x)
		break
		case "ln":
			z=ln(x)
		break
		case "tanh":
			z=tanh(x)
		break
		case "sq":
			z=x*x
		break
		case "sqrt":
			z=sqrt(abs(x))
		break
		case "amoeba":
			z=amoeba(x,1,1)
		break
		case "r":
			z=1/x
		break;
		case "rsq":
			z=1/(x*x)
		break
		case "rsqrt":
			z=1/sqrt(abs(x))
		break
		default:
			z=x
	endswitch
else
	strswitch(fun)
		case "log":
			z=10^x
		break
		case "ln":
			z=exp(x)
		break
		case "tanh":
			z=atanh(x)
		break
		case "sq":
			z=sqrt(abs(x))
		break
		case "sqrt":
			z=x*x
		break
		case "amoeba":
			z=amoeba(x,1,-1)
		break
		case "r":
			z=1/x
		break
		case "rsq":
			z=1/sqrt(abs(x))
		break
		case "rsqrt":
			z=1/(x*x)
		break
		default:
			z=x
	endswitch
endif
return z
end

function amoeba(x,z,d)
Variable x,z,d
if(d>0)
	return (x*x-z*z)/x
else
	return 0.5*(x+sqrt(x*x+4*z*z))
endif
end


//================ Logistic==============================


Proc yBinary(t)  // the first starred y
Variable t=.5
Prompt t,"threshold"
Silent 1
PauseUpdate
Variable /G reg_type=1 // logistic regression
choose_y(yname_set)
Variable k=ysubset[0]
String s=replacesymbols(yname_set[k],"*","")
sprintf s,"*binary(%s).y1",removefileextension(s)
Variable n=DimSize(yval_set,1)
Make /N=(1,n) /O /R yval_set_
Make /N=1 /O /T yname_set
Variable i=0,z
status_set=1
do
	z=yval_set[k][i]
	if(!numtype(z))
		if(z>=t)
			z=1
		else
			z=0
		endif
		yval_set_[0][i]=z
	else
		status_set[i]=2
		yval_set_[0][i]=Nan
	endif	
i+=1
while(i<n)
yname_set[0]=s
duplicate /O yval_set_  yval_set
KillWaves yval_set_
UnstarAllDescriptors()
Make /N=1 /O /I ysubset={0}
end 

function EvSubset_LGR(v,opt) 
Variable opt
Wave /T v
Wave dval_set,yval_set
Wave /I subset_train,subset_valid
Variable n
sort v,v
purge_repeats(v)
//
switch(opt)
	case 0: // training
		n=numpnts(subset_train)
		if(n)
			Make /N=(1,n) /O yv
			yv[0][]=round(yval_set[0][subset_train[q]])
			AutoFill(v,xv,subset_train,0)
			 if(LogisticReg_(xv,yv)>0)
				return Logistic_r2()
			else
				return 1e3
			endif
		endif
	break
	case 1: // validation
		n=numpnts(subset_valid)
		if(n)
			Make /N=(1,n) /O yvb
			yvb[0][]=round(yval_set[0][subset_valid[q]])
			AutoFill(v,xvb,subset_valid,1)
			return LogisticRegB_(xvb,yvb)
		endif
	break
	case -1: // everything
		n=DimSize(dval_set,1)
		if(n)
			Make /N=(1,n) /O yv
			Make /N=(n) /O /I subset_train
			subset_train=p
			yv[0][]=round(yval_set[0][q])
			AutoFill(v,xv,subset_train,0)
			 if(LogisticReg_(xv,yv)>0)
				return Logistic_r2()
			else
				return 1e3
			endif
		endif
	break
endswitch
return 1e3
end


Macro LGR_Genetics(n,l)
Variable n=Ng,l=1e-3
Prompt n,"number of genes"
Prompt l,"regularization, lambda (1e-3 by default)"
Silent 1
PauseUpdate
Variable /G lambda=max(1e-4,l),lev=nan
Variable /G reg_type=1  // logistic regression
Prepare_LGR_genetics()
make_random_vector()
Variable /G Ng=n
Genetic(NT)
hat_eigen_LGR()
end

Proc Prepare_LGR_genetics()
Silent 1
PauseUpdate
avg_std()
Ng=min(Ng,floor(DimSize(dval_set,0)*.6))
printf "%d random descriptors\r\n",Ng
Make /N=1 /O /I ysubset={0}
Make /O subset={0}
Variable /G yshow=0,r2show=nan
Make /N=0 /O /R yp1,yv1,pv1
Make /N=0 /O M_product,M_A,M_B
Make /N=(0,0) /O /R xv,yv,xvb,yvb
	mkucomp()
	mkparam()
//
Prepare_regression()
Prepare_Genetics(5e4)
DoWindow /K Table_coeff
DoWindow /K GeneticsPanel
GeneticsPanel()
end



// opt1: 0 - trainging set, 1 - validation set
//
function AutoFill(v,xv,subset,opt1)
Wave v,xv
Wave /I subset
Variable opt1
Variable /G auto_quad
if(!numpnts(subset))
	return 0
endif
if(!auto_quad)
	return AutoLin(v,xv,subset,opt1)	
else
	return AutoQuad(v,0.9,xv,subset,opt1)
endif
end

function purge_repeats(v)
Wave v
Variable i,n=numpnts(v)
if(n<2)
	return n
endif
Variable flag=0
Variable z=v[0]
Variable m=1
for(i=1;i<n;i+=1) 
	if(v[i]!=z)
		v[m]=v[i]
		z=v[i]
		m+=1
	endif
endfor
DeletePoints m,n-m,v
return m
end
		
function AutoLin(v,xv,subset,opt1)
Wave v,xv
Wave /I subset
Variable opt1
Wave dval_set
Variable j,u
Variable m=numpnts(v)
Variable n=numpnts(subset)
Variable /G reg_type
Variable opt=reg_type
Redimension /N=(m,n) xv
Make /N=(m) /O /T xn
Make /N=(m) /O /B hash
//
for(j=0;j<m;j+=1)
	u=v[j]
	xv[j][]=dval_set[u][subset[q]]
	xn[j]=num2str(u)
endfor
if(!opt1)  // generate hash for parsing - only for training sets
	hash=0
endif
AutoShiftX(xv,opt1)
return AutoTrimX(xv,xn,opt,opt1)
end

function AutoQuad(v,ro,xv,subset,opt1)
Wave v,xv
Wave /I subset
Variable ro,opt1
Wave dval_set
Variable /G reg_type
Variable opt=reg_type
Variable m=numpnts(v)
Variable n=numpnts(subset)
Variable m2=m+m*(m+1)/2
Redimension /N=(m2,n) xv
Make /N=(m2) /O /T xn
Make /N=(m2) /O /B hash
Variable k,i,j,j1,u,u1,z,z1,mm,i1
String s
for(i=0;i<n;i+=1)
	k=m
	i1=subset[i]
	for(j=0;j<m;j+=1)
		u=v[j]
		z=dval_set[u][i1]
		xv[j][i]=z
		if(!i)
			sprintf s,"%s",num2str(u)
			xn[j]=s
		endif
		for(j1=0;j1<=j;j1+=1)
			u1=v[j1]
			z1=dval_set[u1][i1]
			xv[k][i]=z*z1
			if(!i)
				sprintf s,"%s,%s",num2str(u),num2str(u1)
				xn[k]=s
			endif
		k+=1
		endfor
	endfor
endfor
//
if(!opt1)  // generate hash for parsing - only for training sets
	hash=0
endif
AutoShiftX(xv,opt1)
AutoParseX(xv,ro,opt1)
return AutoTrimX(xv,xn,opt,opt1)
end

function AutoTrimX(xv,xn,opt,opt1)
Wave xv
Wave /T xn
Variable opt,opt1
Wave /B hash
Wave xav,xstd
String /G nindex=""
Variable m=DimSize(xv,0)
Variable n=DimSize(xv,1)
Variable j,k,m1=m
for(j=0;j<m1;j+=1)
	m-=hash[j]
endfor
if(!m)
	return 0
endif
Make /N=(m+opt,n) /O /R xv_
if(opt)
	xv_[0][]=1
endif
k=0
for(j=0;j<m1;j+=1)
	if(!hash[j])
		xv_[k+opt][]=xv[j][q]
		nindex+=xn[j]+";"
		k+=1
	endif
endfor
Redimension /N=(m+opt,n) xv
xv=xv_
Redimension /N=0 xv_
//
if(!opt1)
	k=0
	Make /N=(m) /O /R xav_,xstd_
	for(j=0;j<m1;j+=1)
		if(!hash[j])
			xav_[k]=xav[j]
			xstd_[k]=xstd[j]
			k+=1
		endif
	endfor
	Redimension /N=(m) xav,xstd
	xav=xav_
	xstd=xstd_
	Redimension /N=0 xav_,xstd_
endif
return 1
end

function AutoShiftX(xv,opt1)
Wave xv
Variable opt1
Wave /B hash
Wave xav,xstd
Variable m=DimSize(xv,0)
Variable n=DimSize(xv,1)
Variable score=0,i,j,z,z1
for(j=0;j<m;j+=1)
	if(!opt1)
		Make /N=(m) /O xav,xstd
		z=0;z1=0
		for(i=0;i<n;i+=1)
			z+=xv[j][i]
			z1+=xv[j][i]^2
		endfor
		z/=n
		z1=sqrt(abs(z1/n-z*z))
		xav[j]=z
		xstd[j]=z1
		if(z1<1e-6)
			hash[j]=1
			score+=1
		endif
	endif
endfor
xv=(xv-xav[p])/xstd[p]
return score
end

function AutoParseX(xv,ro,opt1)
Wave xv
Variable ro,opt1
Wave /B hash
Variable m=DimSize(xv,0)
Variable n=DimSize(xv,1)
Variable score=0,i,j,j1,z1
for(j=0;j<m;j+=1)
	if(hash[j])
		continue
	endif
	for(j1=j+1;j1<m;j1+=1)
		if(hash[j1])
			continue
		endif
		z1=0
		for(i=0;i<n;i+=1)
			z1+=xv[j][i]*xv[j1][i]
		endfor
		if(sqrt(z1/n)>ro)
			if(!opt1)
				hash[j1]=1
			endif
			score+=1
		endif
	endfor
endfor
return score
end

function /S index2name(s)
String s
Wave /T dname_set
Variable j=index2ind(s,0)
Variable j1=index2ind(s,1)
if(j1<0)	
	s=dname_set[j]
else
	if(j==j1)
		sprintf s,"{%s}^2",dname_set[j]
	else
		sprintf s,"{%s}{%s}",dname_set[j],dname_set[j1]
	endif
	replacesymbols(s,"*","")
endif
return s
end

function index2ind(s,k)
String s
Variable k
Variable n=ItemsInList(s,",")
if(k>=n)
	return -1
else
	return str2num(StringFromList(k,s,","))
endif
end


function accept_best() // set needs to be unshifted
Wave /T dname_set
Wave dval_set
SVAR nindex = nindex
Variable m=ItemsInList(nindex)
Variable n=DimSize(dval_set,1)
Variable j,u1,u2,z,z1,i
String s
Make /N=(m) /O /T dname_set_
Make /N=(m,n) /O /R dval_set_
for(j=0;j<m;j+=1)
	s=StringFromList(j,nindex)
	dname_set_[j]=index2name(s)
	u1=index2ind(s,0)
	u2=index2ind(s,1)
	if(u2<0)
		dval_set_[j][]=dval_set[u1][q]
	else
		dval_set_[j][]=dval_set[u1][q]*dval_set[u2][q]
	endif
endfor
return m
end

function mk_beta0()
Wave beta,xav,xstd,yav,ystd
Variable ny=DimSize(beta,0)
Variable m=DimSize(beta,1)
Make /N=(ny,m+1) /O beta0
Variable i,j,u,z
for(i=0;i<ny;i+=1)
	z=yav[i]
	for(j=0;j<m;j+=1)
		u=ystd[i]*beta[i][j]/xstd[j]
		beta0[i][j+1]=u
		z-=xav[j]*u
	endfor
beta0[i][0]=z
endfor
end


Window bar_test() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(9.75,46.25,463.5,521) p_test_show
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph msize=3
	ModifyGraph grid=2
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=33
	ModifyGraph axOffset(bottom)=1.75
	ModifyGraph axThick=2
	ModifyGraph lblPos(left)=69
	ModifyGraph lblLatPos(bottom)=-9
	Label left "predicted values"
	Label bottom "compounds"
	SetAxis left 0,1.05
	Cursor/P A p_test_show 26
	ShowInfo
	ShowTools
	SetVariable setvar1,pos={240,577},size={90,16},proc=SetYShowProc,title="r2"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={111,577},size={120,16},title=" ",frame=0
	SetVariable setvar5,value= yname_show
	Button button1,pos={502,575},size={25,20},proc=ShiftRight,title=">>"
	Button button2,pos={464,575},size={25,20},proc=ShiftLeft,title="<<"
	SetVariable setvar0,pos={339,577},size={90,16},proc=SetPShowProc,title="show y"
	SetVariable setvar0,value= yshow
EndMacro


Proc SetPShowProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
Variable /G yshow
yshow=CheckCircular(yshow,DimSize(p_test,0))
p_test_show=p_test[yshow][p]
end


Window Table_test() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,42.5,479.25,581) comp_test,smi_test,p_test
	ModifyTable width(Point)=26,alignment(comp_test)=0,width(comp_test)=72,rgb(comp_test)=(0,0,65280)
	ModifyTable alignment(smi_test)=0,width(smi_test)=258,format(p_test)=3,width(p_test)=56
	ModifyTable rgb(p_test)=(65280,0,0),elements(p_test)=(-3,-2)
EndMacro


/// ================== show multiple y's ================



Window multi_y_MMLR_plot() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(78,51,632,682) yp0,yv0
	AppendToGraph/L=new ep0
	ModifyGraph mode(yv0)=3,mode(ep0)=1
	ModifyGraph marker(yp0)=19,marker(yv0)=19
	ModifyGraph lSize(yp0)=2,lSize(yv0)=3,lSize(ep0)=2
	ModifyGraph rgb(yv0)=(0,0,0),rgb(ep0)=(0,39168,0)
	ModifyGraph msize(yv0)=3
	ModifyGraph opaque(yv0)=1
	ModifyGraph grid(left)=2,grid(bottom)=2
	ModifyGraph zero(new)=4
	ModifyGraph mirror(left)=2,mirror(bottom)=2
	ModifyGraph minor(left)=1,minor(bottom)=1
	ModifyGraph fSize=16
	ModifyGraph fStyle=1
	ModifyGraph lblMargin(bottom)=28
	ModifyGraph axOffset(bottom)=1.41667
	ModifyGraph axThick=2
	ModifyGraph notation(left)=1,notation(new)=1
	ModifyGraph zeroThick(new)=1
	ModifyGraph lblPos(left)=69,lblPos(new)=65
	ModifyGraph lblLatPos(bottom)=-10,lblLatPos(new)=-5
	ModifyGraph freePos(new)=2
	ModifyGraph axisEnab(left)={0,0.7}
	ModifyGraph axisEnab(new)={0.75,1}
	Label left "y's"
	Label new "res."
	SetAxis/A/N=1 left
	SetAxis/N=1 bottom 0,35
	SetAxis/A/N=1/E=2 new
	Cursor/P A yp0 7
	ShowInfo
	Legend/C/N=text0/J/A=MC/X=42.39/Y=15.81 "\\Z09\r\\s(yv0) exp \r\\s(yp0) fit "
	ShowTools/A
	SetVariable setvar1,pos={336.00,602.00},size={120.00,19.00},proc=SetYShowProc,title="Pearson ro"
	SetVariable setvar1,format="%.4f",value= r2show
	SetVariable setvar5,pos={206.00,602.00},size={120.00,19.00},title=" ",frame=0
	SetVariable setvar5,value= show_comp
	SetVariable setvar0,pos={97.00,603.00},size={100.00,19.00},proc=SetMultiYShowProc,title="compound"
	SetVariable setvar0,value= show_comp_index
EndMacro


function which_comp_index(q)
Variable q
Wave /B status_set
Variable n=numpnts(status_set)
Variable i,m1=0,m2=0
//
for(i=0;i<=q;i+=1)
	if(status_set[i]==1)	
		m1+=1
	else
		m2+=1
	endif
endfor
if(status_set[q]==1)
	return m1
else
	return -m2
endif
end

function r2vp_to_r2y()
Wave r2vp,ysubset
Wave /T yname_set
Variable m=numpnts(ysubset)
Variable n=numpnts(yname_set)
Make /O /N=(n) /R r2y
r2y=nan
Variable i,j
for(i=0;i<m;i+=1)
	j=ysubset[i]
	r2y[j]=r2vp[i]
endfor
end



Proc SetMultiYShowProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
Silent 1
PauseUpdate
Variable /G r2show=nan
Variable n=DimSize(yv,0)
Make /N=(n) /O /R ep0,yp0,yv0
Variable m=DimSize(dval_set,1)
Variable /G show_comp_index
if(show_comp_index>=m)
	show_comp_index=0
endif
if(show_comp_index<0)
	show_comp_index=m-1
endif
Variable i=which_comp_index(show_comp_index)
Variable i1
String /G show_comp=comp_set[show_comp_index]
yv0=yval_set[ysubset[p]][show_comp_index]
unshift_y()
if(i>0)
	i1=i-1		
	yp0=yp[p][i1]		
	ModifyGraph rgb(yp0)=(65535,0,0)
else
	i1=-i-1
	yp0=ypb[p][i1]
	ModifyGraph rgb(yp0)=(0,0,65535)
endif
ep0=yp0-yv0
r2show=Correlation_1D(yv0,yp0)
end



Proc ClassifyDescriptors()
printf "\r\n%d unhashed, %d starred descriptors\r\n",count_unhashed_des(),count_starred_des()
end

function count_unhashed_des()
Wave /T dname_set
Variable i,m=0,n=numpnts(dname_set)
String s
for(i=0;i<n;i+=1)
	s=dname_set[i]
	if(strsearch(s,"#",0)!=0)
		m+=1
	endif
endfor
return m
end

function count_starred_des()
Wave /T dname_set
Variable i,m=0,n=numpnts(dname_set)
String s
for(i=0;i<n;i+=1)
	s=dname_set[i]
	if(strsearch(s,"*",0)==0)
		m+=1
	endif
endfor
return m
end

//==================================================


Proc Copy_des_wave(s)
String s="natom"
Prompt s,"descriptor"
Variable /G nset
Variable i=find_des_index(s)
if(i<0)
	print "cannot find the descriptor"
else
	s+="_copy"
	Make /O /N=(nset) /R $s
	$s=dval_set[i][p]
endif
end

Proc Hist_des_wave(s)
String s
PauseUpdate
String q=s+"_hist"
Copy_des_wave(s)
s+="_copy"
if(exists(s)==1)
	Make/N=50/O $q
	Histogram/P/C/B=1 $s,$q
	$q*=100
endif
end


Proc Copy_y_wave(s)
String s="EA"
Prompt s,"target (y)"
Variable /G ny
Variable i=find_y_index(s)
if(i<0)
	print "cannot find the target"
else
	s=replacesymbols(s,"*","")
	s+="_copy"
	Make /O /N=(nset) /R $s
	$s=yval_set[i][p]
endif
end

Proc Hist_y_wave(s)
String s
Variable n=50
PauseUpdate
String q=s+"_hist"
Copy_y_wave(s)
s+="_copy"
if(exists(s)==1)
	Make/N=50/O $q
	//Histogram/P/C/B=1 $s,$q
	zhist($s,q)
	$q*=100
endif
end

function zhist(w,s)
Wave w
String s
Wave h=$s
Variable i,j,n=WaveMax(w),m=numpnts(w)
Make /N=(n+1) /O /R $s
h=0
for(i=0;i<m;i+=1)
	j=w[i]
	h[j]+=1
endfor
h/=m
end


function find_des_index(s)
String s
Variable /G ndset
Variable i
String q
Wave /T dname_set
for(i=0;i<ndset;i+=1)
	q=dname_set[i]
	q=replacesymbols(q,"*","")
	q=replacesymbols(q,"#","")
	if(stringmatch(q,s))
		return i
	endif
endfor
return -1
end

function find_y_index(s)
String s
Variable /G ny
Variable i
String q
Wave /T yname_set
for(i=0;i<ny;i+=1)
	q=StringFromList(0,yname_set[i],".")
	q=replacesymbols(q,"*","")
	q=replacesymbols(q,"#","")
	if(stringmatch(q,s))
		return i
	endif
endfor
return -1
end



Proc ExportL_LGR_result(th,opt)
Variable th=-0.2, opt=1
Prompt opt,"0 - export all, 1 - only satisfying the condition"
Prompt th,"log odds threshold for export, + for >, - for <"
Variable /G refnum
String s
Close /A
if(!opt)
	s="logistic_all.csv"
else
	s="logistic_selection.csv"
endif
Open /P=home refnum as s
export_log(th,opt)
Close /A
end

function export_log(th,opt)
Variable th,opt
Variable /G refnum
Wave /T yname_set,comp_set,smi_set
Wave yv1,pv1,status_set
Variable z=(th>0) ? 1 : -1
th*=z
Variable i,n=numpnts(comp_set)
fprintf refnum,"name,smiles,status,%s,log odds",yname_set[0]
for(i=0;i<n;i+=1)
	if(opt && (pv1[i]-th)*z<0)
		continue
	endif
	fprintf refnum,"\r\n%s,%s,%d,%g,%g",comp_set[i],smi_set[i],status_set[i],yv1[i],pv1[i]
endfor
end

//=========================== fit errors and neighbor analysis ==============

Proc ErrorNeighborSort(cutoff)
Variable cutoff = 1.5 
Prompt cutoff,"cutoff error/rms"
deviations_sort(cutoff)
DoWindow /K Table_error_neighbors
if(!V_flag)
	Table_error_neighbors()
endif
end


function deviations_sort(cutoff)
Variable cutoff
Wave /T comp_set, smi_set
Wave ep1,status_set
Variable i,j,n=numpnts(ep1),m=0,eps,k,l
Make /O /N=0 /T comp_set_sorted
String s,q,list=""
for(i=0;i<n;i+=1)
	s = comp_set[i]
	eps = ep1[i]
	if(!numtype(eps) && eps>cutoff && status_set[i]==1)
		j=WhereInSWave(comp_set_sorted,s) 
		if(j<0)
			Make /O /N=(m+1) /T comp_set_sorted,neighbors_set_sorted
			Make /O /N=(m+1) /R ep1_sorted
			ep1_sorted[m]=ep1[i]
			comp_set_sorted[m]=s
			q = list_neighbors(i,1)
			neighbors_set_sorted[m] = q
			for(k=0;k<ItemsInList(q);k+=1)
				s = StringFromList(k,q)
				if(WhichListItem(s,list)<0)
					list+=s+";"
				endif
			endfor
			m+=1
		else
			if(eps>ep1_sorted[j])
				ep1_sorted[j]=eps
			endif
		endif
	endif	
endfor
//
Sort /R ep1_sorted,ep1_sorted,comp_set_sorted,neighbors_set_sorted
list = SortList(list,";",16)
n=ItemsInList(list)
Variable /G refnum
Open /P=home refnum as "neighbor_error_sort.csv"
fprintf refnum,"neighbor,smiles"
for(i=0;i<n;i+=1)
	s = StringFromList(i,list)
	if(strlen(s))
		j = WhereInSWave(comp_set,s) 
		q = smi_set[j]
		print s,q
		fprintf refnum,"\n%s,%s",s,q
	endif
endfor
Close /A
end


function coef_wt_distance(a,b)
Variable a,b
Wave VB,beta1,dval_set,std_set
Wave /T dname_set
Variable i,j,np=numpnts(VB)
Variable x,d=0,d0=0
for(i=0;i<np;i+=1)
	j=VB[i]
	if(strsearch(dname_set[j],".par",0)<0)
		x=dval_set[j][a]-dval_set[j][b]
		d+=(beta1[i]*(x/std_set[j]))^2
		d0+=beta1[i]^2
	endif
endfor
return sqrt(d/d0)
end

// computes the nearest neighbors
function /S list_neighbors(a,opt) // if opt=1 only for status!=1
Variable a,opt
Variable order=5 // how many nearest neighbors
Variable i,m=0,n=numpnts(comp_set)
Wave /T comp_set
Wave status_set
String s,sa=comp_set[a]
Make /O /N=0 /T neighbor_set
for(i=0;i<n;i+=1)
	s = comp_set[i]
	if(opt && status_set[i]==1)
		continue
	endif
	if(!stringmatch(s,sa) && WhereInSWave(neighbor_set,s)<0) 
		Make /O /N=(m+1) /T neighbor_set
		Make /O /N=(m+1) /R neighbor_dist		
		neighbor_set[m] = s
		neighbor_dist[m] = coef_wt_distance(a,i)
		m+=1
	endif
endfor
//
Sort neighbor_dist, neighbor_dist,neighbor_set
n=numpnts(neighbor_set)
s=""
for(i=0;i<min(n,order);i+=1)
	s+=neighbor_set[i]+";"
endfor
return s
end

Window Table_error_neighbors() : Table
	PauseUpdate; Silent 1		// building window...
	Edit/W=(5.25,44.75,660.75,413.75) comp_set_sorted,ep1_sorted,neighbors_set_sorted
	ModifyTable format(Point)=1,width(Point)=24,alignment(comp_set_sorted)=0,rgb(comp_set_sorted)=(1,4,52428)
	ModifyTable alignment(ep1_sorted)=0,format(ep1_sorted)=3,digits(ep1_sorted)=2,width(ep1_sorted)=65
	ModifyTable alignment(neighbors_set_sorted)=0,width(neighbors_set_sorted)=285
EndMacro

Proc ReplaceNansInDes(s,val)
	Variable val=1
	String s = "c.par1"
	Prompt s,"descriptor"
	Prompt val,"value to replace Nan"
	replace_nans_des(s,val)
end

function replace_nans_des(s,val)
	String s
	Variable val
	Wave /T comp_set
	Wave status_set,dval_set
	Variable i, m, n=numpnts(comp_set), q=0
	m = find_des_index(s)
	if (m<0)
		print("No such descriptor")
	else
		for(i=0;i<n;i+=1)
			if(status_set[i]!=1 && numtype(dval_set[m][i]))
				dval_set[m][i] = val
				q+=1
			endif
		endfor
	   printf "\rreplaced %d Nans in %s with %.3f\r",q,s,val
	 endif
end


function find_molecule(s)
String s
Wave /T comp_set
Variable i,n=numpnts(comp_set)
for(i=0;i<n;i+=1)
	if(strsearch(comp_set[i],s,0)>=0)
		print i
		break	
	endif
endfor
end


