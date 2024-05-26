plist := [[0, -120], [0, -1250], [0, -4375], [0, -6250], [1, -3125]];;
selfdeg := 4;;

mono_opr := function(deg, coeff, termnum, clist)
	local cdn, num1, num2;;
	for num1 in [0 .. selfdeg] do
	cdn := coeff * deg^num1;;
	if cdn <> 0 then
		for num2 in [0 .. Length(plist[num1+1]) - 1] do
		if deg + num2 < termnum then
			clist[deg+num2+1] := clist[deg+num2+1] + plist[num1+1][num2+1] * cdn;;
		else
			break;;
		fi;;
		od;;
	fi;;
	od;;
end;;

log_mono_opr := function(deg, logdeg, coeff, termnum, clist)
	local cdn, num1, num2;;
	for num1 in [0 .. selfdeg] do
		if deg <> 0 or num1 >= logdeg then
			cdn := coeff * deg^(num1 - logdeg);;
			for num2 in List([0 .. logdeg - 1], i -> num1 - i) do
				cdn := cdn * num2;;
			od;;
		else
			cdn := 0;;
		fi;;
		if cdn <> 0 then
			for num2 in [0 .. Length(plist[num1+1]) - 1] do
			if deg + num2 < termnum then
				clist[deg+num2+1] := clist[deg+num2+1] + plist[num1+1][num2+1] * cdn;;
			else
				break;;
			fi;;
			od;;
		fi;;
	od;;
end;;

hol_sol := function(inlst, termnum)

local clist, outlst, den, newterm, num, num1;; 
clist := List([1..termnum], num->0);;
outlst := [];;

for num in [0 .. termnum - 1] do
	den := 0;;
	for num1 in [0 .. selfdeg] do
	den := den + plist[num1+1][1] * num^num1;;
	od;;
	
	if den = 0 then
		if num < Length(inlst) then
			Add(outlst, inlst[num+1]);;
			mono_opr(num, inlst[num+1], termnum, clist);;
		else
			Print(1/0);;
		fi;;
	else
		newterm := -clist[num+1] / den;;
		Add(outlst, newterm);;
		mono_opr(num, newterm, termnum, clist);;
	fi;;
od;;

return outlst;
end;;

log_sol := function(inlst, termnum)

local clist, outlst, den, newterm, num, num1, num2;; 
clist := List([1..termnum], num->0);;
outlst := [];;

for num2 in [1 .. Length(inlst)-1] do
	for num in [0 .. Length(inlst[num2+1])-1] do
		if num < termnum then
			log_mono_opr(num, num2, inlst[num2+1][num+1], termnum, clist);;
		else
			break;;
		fi;;
	od;;
od;;

for num in [0 .. termnum - 1] do
	den := 0;;
	for num1 in [0 .. selfdeg] do
		den := den + plist[num1+1][1] * num^num1;;
	od;;
	
	if den = 0 then
		if num < Length(inlst[1]) then
			Add(outlst, inlst[1][num+1]);;
			mono_opr(num, inlst[1][num+1], termnum, clist);;
		else
			Print(1/0);;
		fi;;
	else
		newterm := -clist[num+1] / den;;
		Add(outlst, newterm);;
		mono_opr(num, newterm, termnum, clist);;
	fi;;
od;;

return outlst;
end;;

all_sol := function(termnum)

local initval, sollist, logdeg, num;;
initval := [1];;
sollist := [hol_sol(initval, termnum)];;

for logdeg in [1 .. 3] do
	initval := [[0, 1]];;
	for num in [1 .. Length(sollist)] do
		Add(initval, sollist[Length(sollist) + 1 - num] / Factorial(num));;
	od;;
	Add(sollist, log_sol(initval, termnum));;
od;;
return sollist;;
end;;

all_sol(10);
