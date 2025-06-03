function FindInFile(s)
String s
Variable /G refnum
String /G slast
Variable /G pos_file_last=CurrPos()
String q=""
do
	FReadLine refnum,slast
	if(StrCompCaseIrr(slast,s))
		pos_file_last=CurrPos()
		return 1
	endif
	if(stringmatch(slast,q))
		if(CheckEOF())
			break
		endif
	else
		q=slast
	endif	
while(1)
	FSetPos refnum,pos_file_last
	return 0
end

function check_io_list(slist)
String slist
String /G slast
Variable /G pos_file_last
Variable i,n=ItemsInList(slist)
String s
for(i=0;i<n;i+=1)
		s=StringFromList(i,slist)
		if(StrCompCaseIrr(slast,s))
			pos_file_last=CurrPos()
			return i+1
		endif
	endfor
return 0
end

function FindInFileList(slist)
String slist
Variable /G refnum
String /G slast
Variable /G pos_file_last=CurrPos()
String s,q=""
Variable u
do
	FReadLine refnum,slast
	u=check_io_list(slist)
	if(u)
		return u
	endif
	if(stringmatch(slast,q))
		if(CheckEOF())
			break
		endif
	else
		q=slast
	endif	
while(1)
	FSetPos refnum,pos_file_last
	return 0
end

function FindInFileListG98(slist)
String slist
Variable i,n=ItemsInList(slist)
Variable /G refnum
String /G slast
Variable /G pos_file_last=CurrPos()
String s,q=""
do
	FReadLine refnum,slast
	if(FindEndFileG98(slast))
		break
	endif
	for(i=0;i<n;i+=1)
		s=StringFromList(i,slist)
		if(StrCompCaseIrr(slast,s))
			pos_file_last=CurrPos()
			return i+1
		endif		
	endfor
	if(stringmatch(slast,q))
		if(CheckEOF())
			break
		endif
	else
		q=slast
	endif	
while(1)
	FSetPos refnum,pos_file_last
	return 0
end

function CheckEOF()
Variable /G refnum
FStatus refnum
if(V_filePos>=V_logEOF)
	return 1
endif
return 0
end

function FindInFileFast(s)
String s
String q=""
Variable /G refnum
String /G slast
do
	FReadLine refnum,slast
	if(StrCompCaseIrr(slast,s))
		return 1
	endif
		if(stringmatch(slast,q))
			if(CheckEOF())
				return 0
			endif
	else
		q=slast
	endif	
while(1)
	return 0
end

function StrCompCaseIrr(s1,s)
string s1,s
if(strsearch(UpperStr(s1),UpperStr(s),0)>=0) 
	return 1
endif
return 0
end

function FindInFileG98(s)
String s
Variable /G refnum
String /G slast
Variable /G pos_file_last=CurrPos()
do
	FReadLine refnum,slast
	if(StrCompCaseIrr(slast,s))
		pos_file_last=CurrPos()
		return 1
	endif
	if(FindEndFileG98(slast))
		//printf "\r\nrecord end detected: %s",slast
		break
	endif
while(1)
	FSetPos refnum,pos_file_last
	return 0
end

function FindInFileG98_test(s)
String s
Variable /G refnum
String /G slast
Variable /G pos_file_last=CurrPos()
do
	FReadLine refnum,slast
	print slast
	if(StrCompCaseIrr(slast,s))
		pos_file_last=CurrPos()
		return 1
	endif
	if(FindEndFileG98(slast))
		break
	endif
while(strlen(slast))
	FSetPos refnum,pos_file_last
	return 0
end


function FindEndFileG98(a)
string a
Make /O /T endfile_G98={"@","Entering Link 1"}
Variable i=0,n=numpnts(endfile_G98)
string s
if(CheckEOF())
	return 1
endif
do
	s=endfile_G98[i]
	if(strsearch(a,s,0)>=0)
		return 1
	endif
i+=1
while(i<n)
return 0
end

function CurrPos()
Variable /G refnum
	FStatus refnum
	return V_filePos
end

function SkipLine(n)
Variable n
Variable /G refnum
String /G slast
Variable i=0
do
	FReadLine refnum,slast
	i+=1
while(i<n)
Variable /G pos_file_last=CurrPos()
end

function FeedUntil(s)
String s
Make /N=0 /O /T data0
Variable /G refnum
String /G slast
String G98_break="----"
Variable n=0
Variable /G pos_file_last
do
	FReadLine refnum,slast
	FStatus refnum
	pos_file_last=CurrPos()	
	//print slast
	if(strsearch(slast,s,0)>=0)
		return 0
	endif
	if(strsearch(slast,G98_break,0)>=0)
		return 0
	endif
	Make /O /T /N=(n+1) data0
	data0[n]=slast[0,strlen(slast)-2]
	n+=1
	if(V_filePos==V_logEOF)
		return -1
	endif
while(1)
end

function FeedUntilList(s)
String s
Make /N=0 /O /T data0
Variable /G refnum
String /G slast
String G98_break="----"
Variable u,n=0
Variable /G pos_file_last
s+=";"+G98_break
do
	FReadLine refnum,slast
	FStatus refnum
	pos_file_last=CurrPos()	
	u=check_io_list(s)
	if(u)
		return u
	endif
	Make /O /T /N=(n+1) data0
	data0[n]=slast[0,strlen(slast)-2]
	n+=1
	if(V_filePos==V_logEOF)
		return -1
	endif
while(1)
end




function FileToString()
String s=""
Variable /G refnum
String /G slast
Variable n=0
do
FReadLine refnum,slast
Make /O /N=(n+1) /T stringfile
stringfile[n]=slast
n+=1
FStatus refnum
while(V_filePos<V_logEOF)
return n
end

Macro ReadFileIntoString(str)
String str
Prompt str,"string name"
Variable /G refnum
Open/R refNum
FileToString()
Close refNum
duplicate /O stringfile $str
End

Function/S ReplaceDotsByPs(s)
	String s
	Variable i=0, len=strlen(s)
	do
		if (char2num(s[i]) == 46)
			s[i,i] = "p"	
		endif
		i += 1
	while(i < len)
	return s
end

function FeedWriteUntil(s,rn1,opt)
String s
Variable rn1,opt
Variable /G refnum, pos_file_last
String /G slast
Variable m,i=0,q=-1
do
	FReadLine  refnum,slast
	m=strsearch(slast,"\r",0)
	if(m>0)
		slast=slast[0,m-1]
	endif
	if(strlen(s))
		q=strsearch(slast,s,0)
	endif		
		if( (opt>=0) || (q<0) )
			fprintf rn1,"\r\n%s",slast
		endif
		if(q>=0)
			break
		endif
		if(opt>0)	
			Make /N=(i+1) /T /O data0
			data0[i]=slast[0,m-1]
		endif
		i+=1
		FStatus refnum
	if(V_filePos==V_logEOF)
		pos_file_last=CurrPos()
		return 0
	endif
while(1)
pos_file_last=CurrPos()
return 1
end
