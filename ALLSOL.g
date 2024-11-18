plist := [[0, -120], [0, -1250], [0, -4375], [0, -6250], [1, -3125]];;
plist := [[0, -200, 20480, -392000, 13392640, -104938496, 151121920, 641875968, 43868684288, 11561861120, -2936454250496, 964350181376, 73259609489408, -17250803253248, -999384934252544, 1856734763155456, 9247532739723264, -95633923656122368, 412750136716820480, -970641718175072256, 900671546962477056, 1504184683355701248, -2324490726420774912, -1069041961547071488, 1945555039024054272], [0, -1200, 59360, -1274880, 25339392, -201297920, -629268480, 4305780736, 89209831424, -44288704512, -5325400834048, 1839370076160, 105630123687936, -45903268282368, -1217829923717120, 2676103927824384, 15142843504918528, -164284629375778816, 693981228004540416, -1556784220116877312, 1114856412053241856, 3292623908817076224, -4518095588671094784, -2300213509679480832, 3891110078048108544], [0, 0, 34544, -701952, 94208, 31778816, -2343223296, 6196723712, 93517250560, -128914030592, -3177078325248, -643297181696, 40002520088576, -6729005793280, -290275364700160, 933388735217664, 8917108020740096, -99737472849674240, 415402296201969664, -882925979045986304, 354271442562449408, 2609149888805470208, -3141471846323453952, -1899956092796928000, 2918332558536081408], [0, -6180, -15488, 435456, 372736, -1400832, -349667328, -1430388736, 67549790208, -97055145984, -1551514992640, 9957680349184, -17800894611456, -8681739517952, 242822418530304, 185877594636288, 1926997206892544, -24296286396088320, 101713621662302208, -207351400324136960, -12798315347312640, 873240930872721408, -892557151149490176, -719450040472436736, 972777519512027136], [25, -1180, -27056, 324608, 5351936, -39143424, -330096640, 2802810880, 130547712, -103168344064, 486673481728, 607393939456, -11352101879808, 35167058001920, 7155952386048, -325625093029888, 648501406990336, -1022906591084544, 7236641936637952, -19506985544187904, -7718571626987520, 98393096546418688, -79164837199872000, -106397541196627968, 121597189939003392]];;
plist := [[0, -125, -5550, 16220, 114136, -338569, 164470, 470378, -768200, 812792, -644368, 1340960, -1414944], [0, -850, -23710, 57120, 254100, -952026, 608202, 224532, -1265472, 1790712, -1419568, 3615968, -4009008], [0, -2285, -37791, 77792, 151398, -776517, 680097, -1051014, 159468, 950616, -605540, 3468000, -4166224], [0, -2870, -24118, 38976, -636, 79818, 67050, -909636, 1009968, -365040, 624376, 1368704, -1886592], [25, -1480, -4531, 2236, 1455, 80304, -18597, -38172, -211116, -218600, 547196, 175712, -314432]];;
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
