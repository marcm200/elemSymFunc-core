/computation takes a few seconds on an AMD Ryzen 9 3900X

period(4
repeatiter(inf
describeon(p4a_description.txt

addfss(1,1,0,2
addfss(2,1,0,2,2
addfss(3,1,0,2,2,2
addfss(4,1,0,2,2,2,2

solve(1,B
del(1
soltall(B,100

solve(102,C
del(102
soltall(C,100

resalg(Bareiss,2
restabvall(A,203,100
del(203

/remove some factors from the result to arrive at a small equation
/describing period-4 components
simplifyoff
setstr(400,1*c^2-1*D
setstr(401,1*c^2+1*c-1*D
setstr(402,1*c^2+2*c-1*D+1
setstr(403,1*c^4-2*c^2*D+4*c*D+1*D^2-1*D
setstr(404,1*c^4+2*c^2*D-1*c*D+1*D^2
setstr(405,1*c^4+2*c^3+2*c^2*D+1*c^2+1*c*D+1*D^2-1*D
setstr(407,1*c^8+4*c^7+6*c^6-1*c^5*D+6*c^5-2*c^4*D^2-4*c^4*D+5*c^4-2*c^3*D^2-5*c^3*D+2*c^3+4*c^2*D^2-3*c^2*D+1*c^2-1*c*D^3+3*c*D^2-2*c*D+1*D^4-2*D^3+1*D^2
simplifyon

defeq(500,304,1,0,400,1,0
del(304
defeq(501,500,1,0,401,1,0
del(500
defeq(502,501,1,0,401,1,0
del(501
defeq(503,502,1,0,401,1,0
del(502
defeq(504,503,1,0,402,1,0
del(503
defeq(505,504,1,0,403,1,0
del(504
defeq(506,505,1,0,404,1,0
del(505
defeq(507,506,1,0,405,1,0
del(506
defeq(600,507,1,0,407,1,0
del(507
del(400,499

print(600
saveeq(600,final_Dc.eq

validate(600

describeoff

save(p4a.data

skriptende
.
