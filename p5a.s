/computation takes a few seconds on an AMD Ryzen 9 3900X
period(5
repeatiter(inf
describeon(p5a_description.txt

addfss(1,1,0,2
addfss(2,1,0,2,2
addfss(3,1,0,2,2,2
addfss(4,1,0,2,2,2,2
addfss(5,1,0,2,2,2,2,2

constrbindiv(10,0,4,1,c,1

solve(1,B
del(1
soltall(B,100

solve(102,D
del(102
soltall(D,100

solve(210,C
del(210
soltall(C,100

elmt(400,304,A,8,303
del(304
del(303

resalg(Bareiss,2
restabvall(A,400,100
del(400

/factor result
simplifyoff
setstr(500,4*c-1
setstr(501,2*c^2+2*c-1
setstr(502,1*c^5+4*c^4+6*c^3-1*c^2*E+5*c^2-3*c*E+3*c+1*E^2-2*E+1
setstr(503,48*c^6+56*c^5-89*c^4-64*c^3*E-28*c^3+96*c^2*E+50*c^2-48*c*E-12*c+8*E
setstr(505,16128*c^19+68544*c^18+7168*c^17*E+99840*c^17+18944*c^16*E+69504*c^16+13056*c^15*E-4736*c^15+80384*c^14*E^2-79616*c^14*E-16032*c^14+123776*c^13*E^2-148224*c^13*E+134816*c^13+34816*c^12*E^3-237760*c^12*E^2+202112*c^12*E-48144*c^12-64512*c^11*E^3-73344*c^11*E^2+40064*c^11*E-98624*c^11-149248*c^10*E^3+335424*c^10*E^2-298368*c^10*E+145088*c^10+43776*c^9*E^4+206208*c^9*E^3-316096*c^9*E^2+301312*c^9*E-126584*c^9-6464*c^8*E^4+640*c^8*E^3+151488*c^8*E^2-171616*c^8*E+71492*c^8+27648*c^7*E^5-86272*c^7*E^4-72832*c^7*E^3-5776*c^7*E^2+72368*c^7*E-24276*c^7-69120*c^6*E^5+102784*c^6*E^4+7104*c^6*E^3-35744*c^6*E^2-23408*c^6*E+4744*c^6+71424*c^5*E^5-62848*c^5*E^4+21024*c^5*E^3+19400*c^5*E^2+5088*c^5*E-492*c^5-39040*c^4*E^5+27648*c^4*E^4-10240*c^4*E^3-4640*c^4*E^2-624*c^4*E+21*c^4+11904*c^3*E^5-9376*c^3*E^4+1888*c^3*E^3+536*c^3*E^2+32*c^3*E-1920*c^2*E^5+2176*c^2*E^4-128*c^2*E^3-24*c^2*E^2+128*c*E^5-288*c*E^4+16*E^4
simplifyon

defeq(700,405,1,0,500,1,0
del(405
defeq(701,700,1,0,500,1,0
del(700
defeq(702,701,1,0,500,1,0
del(701
defeq(703,702,1,0,500,1,0
del(702
defeq(704,703,1,0,501,1,0
del(703
defeq(705,704,1,0,502,1,0
del(704
defeq(706,705,1,0,503,1,0
del(705
defeq(707,706,1,0,503,1,0
del(706
defeq(600,707,1,0,505,-1,0
del(707

del(500,599
del(700,799

save(p5a.data

saveeq(600,final_cE.eq
print(600

validate(600

describeoff

skriptende
.
