/computation takes XXX on an AMD Ryzen 9 3900X
period(6
repeatiter(inf
describeon(p6a_description.txt

addfss(1,1,0,2
addfss(2,1,0,2,2
addfss(3,1,0,2,2,2
addfss(4,1,0,2,2,2,2
addfss(5,1,0,2,2,2,2,2
addfss(6,1,0,2,2,2,2,2,2

constrbindiv(10,0,5,1,c,1
constrbindiv(12,0,5,1,c,2

solve(1,B
del(1
soltall(B,100

solve(102,D
del(102
soltall(D,100

solve(210,C
del(210
soltall(C,100

solve(306,E
del(306
soltall(E,100

/factoring
simplifyoff
setstr(500,4*c-1
setstr(501,2*c^2-3*c*A-4*c+1*A+1
simplifyon

defeq(700,412,1,0,500,1,0
del(412
defeq(612,700,1,0,501,1,0
del(700

defeq(700,403,1,0,501,1,0
del(403
defeq(603,700,1,0,501,1,0
del(700

defeq(700,404,1,0,501,1,0
del(404
defeq(604,700,1,0,501,1,0
del(700

defeq(605,405,1,0,501,1,0
del(405
del(500,599

saveeq(612,612.eq
del(612

elmtall(A,8,603,100
del(603
elmtall(A,7,704,100
del(704

loadeq(812,612.eq

/factoring
simplifyoff
setstr(900,1*c+2
setstr(901,4*c-1
setstr(902,1*c^2+3*c-1
setstr(903,2*c^2+2*c-1
setstr(904,2*c^3-2*c^2*A+1*c^2-2*c*A-4*c+1*A+1
simplifyon

defeq(1100,805,1,0,900,1,0
del(805
defeq(1101,1100,1,0,901,1,0
del(1100
defeq(1102,1101,1,0,902,1,0
del(1101
defeq(1103,1102,1,0,903,1,0
del(1102
defeq(1104,1103,1,0,903,1,0
del(1103
defeq(1105,1104,1,0,903,1,0
del(1104
defeq(1005,1105,1,0,904,1,0
del(1105
del(900,999

resalg(Bareiss,2
restabvall(A,812,300
del(812

/factor result
simplifyoff
setstr(1400,+1*c^3++2*c^2++1*c+-1*F++1
setstr(1401,+1*c^3++3*c^2++3*c+-1*F++1
setstr(1402,+1*c^6++4*c^5++6*c^4++2*c^3*F++6*c^3++3*c^2*F++5*c^2+-2*c*F++2*c++1*F^2+-2*F++1
setstr(1403,+1*c^9++6*c^8++15*c^7++1*c^6*F++22*c^6++3*c^5*F++23*c^5++2*c^4*F++18*c^4+-1*c^3*F^2+-1*c^3*F++10*c^3+-1*c^2*F^2+-4*c^2*F++5*c^2++3*c*F^2+-6*c*F++3*c+-1*F^3++3*F^2+-3*F++1
setstr(1404,+48*c^13++288*c^12++696*c^11++952*c^10++943*c^9++620*c^8+-128*c^7*F^2++64*c^7*F+-10*c^7+-320*c^6*F^2++128*c^6*F+-320*c^6+-96*c^5*F+-167*c^5++352*c^4*F^2+-224*c^4*F++48*c^4+-32*c^3*F^2++80*c^3*F++46*c^3+-144*c^2*F^2++96*c^2*F+-12*c^2++64*c*F^2+-56*c*F+-8*F^2++8*F
simplifyon

defeq(1500,1305,1,0,1400,1,0
del(1305
defeq(1501,1500,1,0,1401,1,0
del(1500
defeq(1502,1501,1,0,1402,1,0
del(1501
defeq(1503,1502,1,0,1403,1,0
del(1502
defeq(1700,1503,1,0,1404,-1,0
del(1503
del(1400,1499

save(p6a.data
saveeq(1700,final_Fc.eq
print(1700
validate(1700

describeoff

skriptende
.
