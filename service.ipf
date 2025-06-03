//=================== service ======================== //

Menu "Service"
	"RemoveTraceCsr"
	"ScaleTraceCsr"
	"NormTraceCsr"
	"IdentifyWaveCursor"
	"Integrate1stDerivativeAB"
	"IntergateAsExperimentAB"
end

Proc Integrate1stDerivativeAB()
Variable xa=xcsr(A)
Variable xb=xcsr(B)
String s=CsrWave(A)
String ss="int_"+s
duplicate /O $s $ss
Integrate $ss
Variable norm=area($ss,-Inf,Inf)
$ss/=norm
print "area to total, % = ",100*area($ss,xa,xb)
end

Proc IntergateAsExperimentAB()
Variable a,b,n
string s=CsrWave(A)
variable norm
a=pcsr(A)
b=pcsr(B)
n=b-a
Make /O /N=(n+1) /R yint
yint=$s[a+p]
variable x1=pnt2x($s,a)
variable x2=pnt2x($s,b)
SetScale/P x x1,deltax($s),"G", yint
variable y1=yint[0]
variable y2=yint[n]
yint-=y1+(y2-y1)*p/(b-a) // subtract linear drift from the derivative
duplicate /O $s dummy
Integrate /T dummy  // first integration overall
norm=area(dummy,-Inf,Inf)
duplicate /O yint dummy 
Integrate /T dummy  // first integration of the highlighted area
dummy/=norm
norm=100*area(dummy,-Inf,Inf) // second integration
print "integral, % of the total =",norm
doWindow /K integration_plot
integration_plot()
EndMacro


Proc RemoveTraceCsr()
string s=CsrWave(A)
remove $s
end

Proc ScaleTraceCsr(a)
Variable a=1
Prompt a,"scaling factor"
string s=CsrWave(A)
$s*=a
end

Proc NormTraceCsr()
Variable a=pcsr(A)
string s=CsrWave(A)
Variable norm=$s[a]
$s/=norm
end

Proc IdentifyWaveCursor()
print CsrWave(A)
end

function TextWaveToFile(w,rn)
Wave /T w
Variable rn
Variable i=0,n=numpnts(w)
do
	fprintf rn,"\r\n%s",w[i]
i+=1
while(i<n)
end

function TextWaveCStyleOut(w)
Wave /T w
Variable rn
Variable i=0,n=numpnts(w)
Variable count=0
do
	printf "\"%s\",",w[i]
	if(count==20)
		printf "\r\n"
		count=0
	endif
	count+=1
i+=1
while(i<n)
printf "\r\n"
end

function WaveCStyleOut(w)
Wave  w
Variable rn
Variable i=0,n=numpnts(w)
Variable count=0
do
	printf "%g,",w[i]
	if(count==20)
		printf "\r\n"
		count=0
	endif
	count+=1
i+=1
while(i<n)
printf "\r\n"
end

function /S replacespacesbyunderbars(s)
String s
return replacesymbols(s," ","_")
end


function /S repeatstring(s,n)
String s
Variable n
Variable i
String q=""
for(i=0;i<n;i+=1)
	q+=s
endfor
return q
end


function /S replacesymbols(s,a,b)
string s,a,b
Variable i=0,n=strlen(s)
Variable k,m=strlen(a)
string c,s1=""
do
	k=strsearch(s,a,i)
	if(k<0)
		break
	endif
	s1+=s[i,k-1]+b
	i=k+m
while(i<n)
return s1+s[i,n-1]
end

function /S removefilepath_(s,symb)
string s,symb
SVAR s_path = root:s_path
Variable i=0,m=0
do
	i=strsearch(s,symb,m)
	if(i<0)
		break
	else
	m=i+1
	endif
while(1)
s_path=s[0,m-1]
	return  s[m,strlen(s)-1]
end



function /S removefilepath(s)
string s
print s
return removefilepath_(s,":")
end

function /S removefileextension(s)
string s
Variable i=strsearch(s,".",0)
	if(i<0)
		return s
	else
		return s[0,i-1]
	endif
end

function /S removefilelastextension(s)
string s
Variable i=0,j=strlen(s)
do
	i=strsearch(s,".",i)
	if(i<0)
		break
	else
		j=i
		i+=1
	endif
while(1)
return s[0,j-1]
end

function /S removefilepathext(s)
string s
s=removefilepath(s)
s=removefileextension(s)
return s
end


function /S replacepoint(s)
string s
string s1
Variable i=strsearch(s,".",0)
if(i<0)
	return s
else
	s[i,i]="p"
	return s
endif
end



function notzero(x)
variable x
if(x)
	return sign(x)
else
	return 0
endif
end

function /S trimCR(s)
String s
Variable n=strsearch(s,"\r",0)
if(n<0)
	return s
else
	return s[0,n-1]
endif
end

function /S path2script()
SVAR s_path = root:s_path
Variable n=strlen(s_path)
String s=s_path[2,n-1]
s="C:\\"+replacesymbols(s,":","\\")
return s
end

function CheckLimits(x,a,b)
Variable x,a,b
Variable a1=min(a,b)
Variable b1=max(a,b)
if(x<a1)
	return a1
endif
if (x>b1)
	return b1
endif
return x
end

function CheckCircular(i,n)
variable i,n
if(i<0)
	return n-1
endif
if(i>=n)
	return 0
endif
return i
end

function CheckCircularBelong(w,x) // w must be sorted
Variable x
Wave w
Variable i=CheckCircularBelongIndex(w,x)
return w[i]
end

function CheckCircularBelongIndex(w,x) // w must be sorted
Variable x
Wave w
Variable i,n=numpnts(w)
for(i=0;i<n-1;i+=1)
	if((x>w[i]) && (x<=w[i+1]))
		return i+1
	endif
endfor
if(x<w[0])
	return n-1
endif	
if(x>w[n+1])
	return 0
endif
return 0
end

function CheckRange(x,a,b)
Variable x,a,b
Variable a1=min(a,b)
Variable b1=max(a,b)
if((x>=a1) && (x<=b1))
	return 1
else
	return 0
endif
end

function /S ChangeStringItem(n,a,s,sep)
Variable n
String a,s,sep
Variable i=0
String sout="",u
Variable nz=ItemsInList(s,sep)
if(n>nz)
	return ""
endif
do
	u=StringFromList(i,s,sep)
	if(i!=n)
		sout+=u
	else
		sout+=a
	endif
	if(i!=nz-1)
		sout+=sep
	endif
i+=1
while(i<nz)
return sout
end

function IsInWave(w,x)
Variable x
Wave W
Variable i,n=numpnts(w)
if(n)
	for(i=0;i<n;i+=1)
		if(x==w[i])
		return 0
	endif
	endfor
endif
return 1
end

function CountInSWave(w,s)
Wave /T w
String s
Variable i,m=0
Variable n=numpnts(w)
for(i=0;i<n;i+=1)
	if(strsearch(w[i],s,0)>=0)
		m+=1
	endif
endfor
return m
end


function WhereInWave(w,x)
Variable x
Wave W
Variable i,n=numpnts(w)
if(!n)
	return -1
endif
for(i=0;i<n;i+=1)
	if(x==w[i])
	return i
endif
endfor
return -1
end

function WhereInSWave(w,s)
String s
Wave /T W
String q
Variable i,n=numpnts(w)
if(!n)
	return -1
endif
s=replacesymbols(s,"*","")
for(i=0;i<n;i+=1)
	q=replacesymbols(w[i],"*","")
	if(stringmatch(s,q))
	return i
endif
endfor
return -1
end

function WhereInSWaveLoose(w,s)
String s
Wave /T W
s=LowerStr(s)
Variable i,n=numpnts(w)
String q
for(i=0;i<n;i+=1)
	q=LowerStr(w[i])
	if((strlen(q)) && (strsearch(q,s,0)>=0))
		return i
	endif
endfor
return -1
end

function CompareSWavesHash(w1,w2,result)
Wave /B result
Wave /T w1,w2
Variable n=numpnts(w1)
Variable m=numpnts(w2)
Variable i,u
if(m<n)
	return -1
endif
Make /N=(n) /O /B result
m=0
for(i=0;i<n;i+=1)
	u=!stringmatch(w1[i],w2[i])
	result[i]=u
	if(u)
		m+=1
	endif
endfor
return m
end

function CompareSWaves(w1,w2)
Wave /T w1,w2
Variable n=numpnts(w1)
Variable m=numpnts(w2)
Variable i,u
if(m<n)
	return -1
endif
Make /N=(n) /O /B result
m=0
for(i=0;i<n;i+=1)
	u=!stringmatch(w1[i],w2[i])
	if(u)
		m+=1
	endif
endfor
return m
end



function /S AnalyzeOrderedList(s)
String s
String s1="",q
s=replacesymbols(s,",",";")
Variable n=ItemsInList(s,";")
Variable i=0,m=0,j,a,b
do
	q=StringFromList(i,s)
	if(strsearch(q,"-",0)>=0)
		sscanf q,"%d-%d",a,b
		j=min(a,b)
		do
			s1+=num2str(j)+";"
		j+=1
		while(j<=max(a,b))
		m=max(m,j)
	else
		sscanf q,"%d",a
		s1+=num2str(a)+";"
		m=max(m,a)
	endif
i+=1	
while(i<n)
i=0
s=""
do
	q=num2str(i)
	if(FindListItem(q,s1,";",0)>=0)
		s+=q+";"
	endif
i+=1
while(i<=m)
return s
end

function GetFromList(i,s)
string s
Variable i
string q=StringFromList(i,s)
return str2num(q)
end

function FindInList(s,list,symb)
String s,list,symb
Variable i,n=ItemsInList(list,symb)
String q
for(i=0;i<n;i+=1)
	q=StringFromList(i,list,symb)
	if(strsearch(q,s,0)>=0)
		return i
	endif
endfor
return -1
end

function FindInListExact(s,list,symb)
String s,list,symb
Variable i,n=ItemsInList(list,symb)
String q
for(i=0;i<n;i+=1)
	q=StringFromList(i,list,symb)
	if(stringmatch(q,s))
		return i
	endif
endfor
return -1
end

function FindInListCast(s,list,symb)
String s,list,symb
Variable i,n=ItemsInList(list,symb)
String q
for(i=0;i<n;i+=1)
	q=StringFromList(i,list,symb)
	if(stringmatch(s,q))
		return i
	endif
endfor
return -1
end

function /S StringFindInList(s,list,symb)
String s,list,symb
Variable i,n=ItemsInList(list,symb)
String q
for(i=0;i<n;i+=1)
	q=StringFromList(i,list,symb)
	if(strsearch(q,s,0)>=0)
		return q
	endif
endfor
return ""
end

function /S StringFindInListExact(s,list,symb)
String s,list,symb
Variable i,n=ItemsInList(list,symb)
String q
for(i=0;i<n;i+=1)
	q=StringFromList(i,list,symb)
	if(stringmatch(s,q))
		return q
	endif
endfor
return ""
end

function similarity(w1,w2,a,b)
Wave w1,w2
Variable a,b
Variable i,x,y
Variable sxy=0,sxx=0
Variable n=numpnts(w1)
Variable /G alpha
for(i=a;i<=b;i+=1)
	x=w1[i]
	y=w2[i]
	if( !(numtype(x)) && !(numtype(y)) )
		sxy+=x*y
		sxx+=x*x
	endif
endfor
return sxy/sxx
end

function Correlation_1D(w1,w2)
Wave w1,w2
Variable i,x,y
Variable sxy=0,sx=0,sy=0,sxx=0,syy=0,m=0
Variable n=numpnts(w1)
Variable /G alpha
for(i=0;i<n;i+=1)
	x=w1[i]
	y=w2[i]
	if( !(numtype(x)) && !(numtype(y)) )
		sx+=x; sy+=y
		sxy+=x*y;sxx+=x*x;syy+=y*y
		m+=1
	endif
endfor
sx/=m
sy/=m
sxx=sxx/m-sx*sx
syy=syy/m-sy*sy
sxy=sxy/m-sx*sy
return sxy/sqrt(sxx*syy)
end

function Correlation_2D(w,a,b,opt)
Wave w
Variable a,b,opt
Variable i,x,y
Variable sxy=0,sx=0,sy=0,sxx=0,syy=0,m=0
Variable n=DimSize(w,opt)
for(i=0;i<n;i+=1)
	if(opt)
		x=w[a][i]; y=w[b][i]
	else
		x=w[i][a]; y=w[i][b]
	endif
	if( !(numtype(x)) && !(numtype(y)) )
		sx+=x; sy+=y
		sxy+=x*y;sxx+=x*x;syy+=y*y
		m+=1
	endif
endfor
sx/=m
sy/=m
sxx=sxx/m-sx*sx
syy=syy/m-sy*sy
sxy=sxy/m-sx*sy
return sxy/sqrt(sxx*syy)
end



// Interquartile range

function IQR_1D(w)
Wave w
duplicate /O w dummy
sort dummy,dummy
Variable n=numpnts(dummy)
Variable q1=dummy[.25*n]
Variable q2=dummy[.5*n]
Variable q3=dummy[.75*n]
return (q3-q1)
end

function IQR_2D(w,a,opt)
Wave w
Variable a,opt
Variable n=DimSize(w,opt)
Make /N=(n) /R /O dummy
if(opt)
	dummy=w[a][p]
else
	dummy=w[p][a]
endif
sort dummy,dummy
Variable q1=dummy[.25*n]
Variable q2=dummy[.5*n]
Variable q3=dummy[.75*n]
return (q3-q1)
end


Proc CopyWaveList(list,s,opt)
String s,list
Variable opt
String a,q
Variable i=0,n=ItemsInList(list)
Silent 1
PauseUpdate
do
	a=StringFromList(i,list)
	q=a+"_"+s
	//print a,q
	if(opt>0)
		if(exists(a)==1)
			duplicate /O $a $q
		endif
	else
		if(exists(q)==1)
			duplicate /O $q $a
		endif
	endif
i+=1
while(i<n)
end

function compare_waves(w1,w2)
Wave w1,w2
Variable i, n1=numpnts(w1), n2=numpnts(w2)
if(n1!=n2)
	return 1
endif
for(i=0;i<n1;i+=1)
	if(w1[i]!=w2[i])
		return 1
	endif
endfor
return 0
end

function WriteLog(name,list)
String name,list
Variable i,n=ItemsInList(list)
Variable /G refnum
if(n)
	Open /P=data refnum as "log_"+name+".txt"
	for(i=0;i<n;i+=1)
		fprintf refnum,"%s\r\n",StringFromList(i,list)
	endfor
	Close refnum
endif
end


function copyval2waves(list,prefix,j)
String list,prefix
Variable j
Variable i=0, n=ItemsInList(list)
String q,s
String fldrSav= GetDataFolder(1)
for(i=0;i<n;i+=1)
	s=StringFromList(i,list)
	sprintf q,"%s_%s",prefix,s
	Make /O /R /N=(j+1) $q
	Wave w=$q
	sprintf q,"%s%s",fldrSav,s
	NVAR x = $q
	w[j]=x
endfor
end

function /S copylist(list,prefix)
String list,prefix
Variable i,n=ItemsInList(list)
String s,q=""
for(i=0;i<n;i+=1)
	s=StringFromList(i,list)
	q+=prefix+"_"+s+";"
	sprintf s,"Variable /G %s=Nan",s
	Execute s
endfor
return q
end

Proc MakeFolder(folder)
String folder
// print "\r\n>> creates folder ",folder
Silent 1
PauseUpdate
Variable i=0,m=ItemsInList(folder,":")
String s,q
do
	q=StringFromList(i,folder,":")
	if(i)
		s+=":"+q
	else
		s=q
	endif
	if(strsearch(q,"root",0)<0)
		NewDataFolder /O $s
	endif
i+=1
while(i<m)
end

Proc copy2folder(list,folder)
String list,folder
Silent 1
PauseUpdate
SetDataFolder root:
folder=replacesymbols(folder,"::",":")
MakeFolder(folder)
Variable i=0,n=ItemsInList(list)
String s,q
do
	s=StringFromList(i,list)
	if(exists(s)==1)
		q=folder+":"+s
		q=replacesymbols(q,"::",":")
		//print "s=",s,"q=",q
		duplicate /O $s $q
	endif
i+=1
while(i<n)
end

Proc copyfromfolder(list,folder)
String list,folder
Silent 1
PauseUpdate
//printf "\r\ncopying from %s to root",folder
Variable i=0,n=ItemsInList(list)
String s,q
SetDataFolder folder
do
	s=StringFromList(i,list)	
	q=folder+":"+s
	s="root:"+s
	// printf "\r\n>> %s -> %s",q,s
	if(exists(q)==1)	
		duplicate /O $q $s
	endif
i+=1
while(i<n)
SetDataFolder root:
end

Proc copy2folderprefix(list,prefix,folder)
String list,folder,prefix
Variable i=0,n=ItemsInList(list)
String s,q
SetDataFolder root:
Execute "NewDataFolder /O "+folder
Silent 1
PauseUpdate
do
	sprintf s,"%s_%s",prefix,StringFromList(i,list)
	if(exists(s)==1)
		sprintf q,"%s:%s",folder,s
		duplicate /O $s $q
	endif
i+=1
while(i<n)
end


Macro Window2Folder()
SetDataFolder root:
NewDataFolder /O root:last_folder
String list=WaveList("*",";","WIN:")
Variable i,w
copy2folder(list,"root:last_folder")
end

Macro WindowSave()
WindowSave_()
end

Macro CsrWaveSave()
WaveSave_()
end

function WindowSave_()
String list=TraceNameList("", ";",1)
Variable n=ItemsInList(list)
Variable i=0
String s,q="",full
for(i=0;i<n;i+=1)
	s=StringFromList(i,list)
	full=GetWavesDataFolder(TraceNameToWaveRef("",s),2)
	q+=full+";"
	wave w=XWaveRefFromTrace("",s)
	if(waveexists(w))
		q+=GetWavesDataFolder(w,2)+";"
		print s,full," vs ",XWaveName("",s)
	else 
		print s,full
	endif
endfor
q=RemoveListDuplicates(q,",")
Execute "Save /T "+q
end

function WaveSave_()
Wave wy=CsrWaveRef(A)
Wave wx=CsrXWaveRef(A)
Save /T wy,wx
end

Function/t RemoveListDuplicates(s,symb)
    String s,symb
    String u,q = ""
    variable i,n=itemsinlist(s)
    for(i = 0 ; i < n ; i+=1)
	u=stringfromlist(i,s)
        if(whichlistitem(u, q) <0 )
            q += u+symb
        endif
    endfor
    return q
end



function VariableFromList(i,s)
Variable i
String s
return str2num(StringFromList(i,s))
end

function digits_filter(s)
String s
String digits="0;1;2;3;4;5;6;7;8;9"
s=UpperStr(s)
String q=""
Variable i,u,n=strlen(s)
for(i=0;i<n;i+=1)
	u=WhichListItem(s[i],digits)
	if(u>=0)
		q+=s[i]
	endif
endfor
if(!strlen(q))
	return Nan
else
	return str2num(q)
endif
end

function /S nodigits_filter(s)
String s
String digits="0;1;2;3;4;5;6;7;8;9"
String q=""
Variable i,u,n=strlen(s)
for(i=0;i<n;i+=1)
	u=WhichListItem(s[i],digits)
	if(u<0)
		q+=s[i]
	endif
endfor
return q
end
 

function /S num2abc(i)
Variable i
String abc="A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z"
Variable n=ItemsInList(abc,",")
i=mod(i,n)
return StringFromList(i,abc,",")
end

function /S num2abcd(i)
Variable i
String abcd="0,1,2,3,4,5,6,7,8,9,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z"
Variable n=ItemsInList(abcd,",")
i=mod(i,n)
return StringFromList(i,abcd,",")
end


function abc2num(s)
String s
s=UpperStr(s[0])
String abc="A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z"
return WhichListItem(s,abc,",")
end

function /S salign(s,n,opt)  // opt=0  left opt=1 right alignmentpri
String s
Variable n,opt
Variable m=strlen(s)
if(m>=n)
	return s[0,n-1]
endif
Variable j
for(j=0;j<n-m;j+=1)
	if(opt)
		s=" "+s
	else
	 	s+=" "
	endif
endfor
return s
end

function appendwave(w,w1)
Wave w,w1
variable n=numpnts(w)
Variable m=n+numpnts(w1)
Redimension /N=(m) w
w[n,m-1]=w1[p-n]
end

function appendval(w,x)
Wave w
Variable x
variable n=numpnts(w)
Redimension /N=(n+1) w
w[n]=x
end

function appendSwave(w,w1)
Wave /T w,w1
variable n=numpnts(w)
Variable m=n+numpnts(w1)
Redimension /N=(m) w
w[n,m-1]=w1[p-n]
end

function appendlist(w,s,code)
Wave /T w
String s,code
variable n=numpnts(w)
variable m=n+ItemsInList(s)
Redimension /N=(m) w
w[n,m-1]=StringFromList(p-n,s)+code
end

function  ButtonIsClicked_()
	GetMouse
	if (V_flag & 2)  // right click
		return 1
	else
		return 0
	endif
end
