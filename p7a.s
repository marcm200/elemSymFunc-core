/computation takes 20 min on an AMD Ryzen 9 3900X
period(7
repeatiter(inf
describeon(p7a_description.txt

addfss(1,1,0,2
addfss(2,1,0,2,2
addfss(3,1,0,2,2,2
addfss(4,1,0,2,2,2,2
addfss(7,1,0,2,2,2,2,2,2,2

constrbindiv(12,0,6,1,c,2
constrbindiv(13,0,6,1,c,3

solve(1,B
del(1
soltall(B,100

solve(102,D
del(102
soltall(D,100

solve(203,F
del(203
soltall(F,100

solve(312,E
del(312
soltall(E,100

elmtall(C,2,413,100
solve(507,C
del(507
soltall(C,100

/factoring
simplifyoff
setstr(700,+4*c+-1
setstr(701,+2*c^2+-3*c*A+-4*c++1*A++1
simplifyon

defeq(900,513,1,0,700,1,0
del(513
defeq(901,900,1,0,700,1,0
del(900
defeq(902,901,1,0,700,1,0
del(901
defeq(903,902,1,0,701,1,0
del(902
defeq(1013,903,1,0,701,1,0
del(903

defeq(900,604,1,0,700,1,0
del(604
defeq(901,900,1,0,700,1,0
del(900
defeq(902,901,1,0,700,1,0
del(901
defeq(903,902,1,0,700,1,0
del(902
defeq(904,903,1,0,701,1,0
del(903
defeq(905,904,1,0,701,1,0
del(904
defeq(1004,905,1,0,701,1,0
del(905

del(700,799

resalg(Bareiss,2
restabvall(A,1013,100
del(1013

simplifyoff
setstr(1200,+1*c++2
setstr(1201,+1*c^7++6*c^6++15*c^5++21*c^4+-1*c^3*G++19*c^3+-4*c^2*G++12*c^2+-5*c*G++5*c++1*G^2+-2*G++1 
setstr(1202,+3*c^7++14*c^6++30*c^5++44*c^4++40*c^3++16*c^2+-8*c*G++8*c++8*G^2+-16*G++8 
setstr(1203,+1*c^21++12*c^20++66*c^19++227*c^18++1*c^17*G++565*c^17++10*c^16*G++1107*c^16++43*c^15*G++1785*c^15++3*c^14*G^2++107*c^14*G++2432*c^14++19*c^13*G^2++176*c^13*G++2863*c^13++40*c^12*G^2++196*c^12*G++2957*c^12++17*c^11*G^2++108*c^11*G++2708*c^11++2*c^10*G^3+-64*c^10*G^2+-93*c^10*G++2228*c^10++15*c^9*G^3+-148*c^9*G^2+-323*c^9*G++1656*c^9++46*c^8*G^3+-183*c^8*G^2+-491*c^8*G++1116*c^8++3*c^7*G^4++77*c^7*G^3+-122*c^7*G^2+-534*c^7*G++688*c^7++7*c^6*G^4++82*c^6*G^3+-3*c^6*G^2+-452*c^6*G++382*c^6+-6*c^5*G^4++59*c^5*G^3++81*c^5*G^2+-327*c^5*G++193*c^5+-18*c^4*G^4++14*c^4*G^3++114*c^4*G^2+-198*c^4*G++88*c^4++1*c^3*G^5+-3*c^3*G^4+-31*c^3*G^3++101*c^3*G^2+-102*c^3*G++34*c^3++2*c^2*G^5++5*c^2*G^4+-40*c^2*G^3++70*c^2*G^2+-50*c^2*G++13*c^2+-3*c*G^5++15*c*G^4+-30*c*G^3++30*c*G^2+-15*c*G++3*c++1*G^6+-6*G^5++15*G^4+-20*G^3++15*G^2+-6*G++1 
setstr(1204,+1*c^21++14*c^20++91*c^19++368*c^18++1*c^17*G++1049*c^17++12*c^16*G++2266*c^16++65*c^15*G++3889*c^15++3*c^14*G^2++213*c^14*G++5474*c^14++23*c^13*G^2++480*c^13*G++6461*c^13++68*c^12*G^2++803*c^12*G++6512*c^12++82*c^11*G^2++1043*c^11*G++5705*c^11++2*c^10*G^3+-39*c^10*G^2++1067*c^10*G++4420*c^10++15*c^9*G^3+-306*c^9*G^2++847*c^9*G++3072*c^9++50*c^8*G^3+-575*c^8*G^2++491*c^8*G++1930*c^8++3*c^7*G^4++101*c^7*G^3+-657*c^7*G^2++150*c^7*G++1091*c^7++9*c^6*G^4++136*c^6*G^3+-480*c^6*G^2+-84*c^6*G++547*c^6++3*c^5*G^4++117*c^5*G^3+-177*c^5*G^2+-189*c^5*G++246*c^5+-11*c^4*G^4++44*c^4*G^3++54*c^4*G^2+-196*c^4*G++109*c^4++1*c^3*G^5+-5*c^3*G^4+-43*c^3*G^3++149*c^3*G^2+-154*c^3*G++52*c^3++23*c^2*G^4+-92*c^2*G^3++138*c^2*G^2+-92*c^2*G++23*c^2+-7*c*G^5++35*c*G^4+-70*c*G^3++70*c*G^2+-35*c*G++7*c++1*G^6+-6*G^5++15*G^4+-20*G^3++15*G^2+-6*G++1 
setstr(1205,+96832*c^22++635968*c^21++1830784*c^20++3292144*c^19++3912377*c^18++2068704*c^17+-342144*c^16*G+-1001992*c^16++404352*c^15*G^2+-1487808*c^15*G+-1783788*c^15++793152*c^14*G^2+-1932336*c^14*G+-494368*c^14++176688*c^13*G^2+-480672*c^13*G++1253040*c^13++14400*c^12*G^2++1461408*c^12*G++424932*c^12+-1221408*c^11*G^2++2534880*c^11*G+-1863328*c^11++156704*c^10*G^2+-1600368*c^10*G++733280*c^10+-839808*c^9*G^3++2739504*c^9*G^2+-2767520*c^9*G++1127312*c^9++419904*c^8*G^4++559872*c^8*G^3+-3092928*c^8*G^2++3600448*c^8*G+-1327232*c^8+-1119744*c^7*G^4++1866240*c^7*G^3+-524480*c^7*G^2+-712032*c^7*G++359840*c^7++1306368*c^6*G^4+-3483648*c^6*G^3++3616864*c^6*G^2+-1767232*c^6*G++364320*c^6+-870912*c^5*G^4++2757888*c^5*G^3+-3385984*c^5*G^2++1932160*c^5*G+-438016*c^5++362880*c^4*G^4+-1257984*c^4*G^3++1662592*c^4*G^2+-997248*c^4*G++230016*c^4+-96768*c^3*G^4++354816*c^3*G^3+-490880*c^3*G^2++304128*c^3*G+-71296*c^3++16128*c^2*G^4+-61440*c^2*G^3++87872*c^2*G^2+-55936*c^2*G++13376*c^2+-1536*c*G^4++6016*c*G^3+-8832*c*G^2++5760*c*G+-1408*c++64*G^4+-256*G^3++384*G^2+-256*G++64 
simplifyon

defeq(1400,1104,1,0,1200,1,0,24
del(1104
defeq(1401,1400,1,0,1201,1,0
del(1400
defeq(1402,1401,1,0,1202,1,0,6
del(1401
defeq(1403,1402,1,0,1203,1,0,1
del(1402
defeq(1404,1403,1,0,1204,1,0,1
del(1403
defeq(1600,1404,1,0,1205,1,0
del(1404
del(1200,1299

save(p7a.data
print(1600
saveeq(1600,final_cG.eq

validate(1600

describeoff

skriptende
.
