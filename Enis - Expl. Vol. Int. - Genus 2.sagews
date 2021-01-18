︠ca180088-915d-4a7a-b321-3c22d3d50f09ss︠
R.<x> = QQ[]
X = HyperellipticCurve(R([5, 0, -4, 10, -8, 0, 1])) # LMFDB 3950.b.39500.1
f = X.hyperelliptic_polynomials()[0]
R, S = X(1,2), X(1,-2)
K = Qp(5,12)
XK = X.change_ring(K)
L.<a> = K.ext(x^4-5) # ramified extension
XL = X.change_ring(L)
P1 = XL.lift_x(a,all)[1] # ref pt in e_1
P2 = XL.lift_x(a,all)[0] # ref pt in e_2
P3 = XL.lift_x(a+2,all)[1] # ref pt in e_3
P4 = XL.lift_x(a+2,all)[0] # ref pt in e_4
P5 = XL.lift_x(a+3,all)[0] # ref pt in e_5
P6 = XL.lift_x(a+3,all)[1] # ref pt in e_6
def Log(z): return L(z).log(p_branch = 0, change_frac = True) # Iwasawa branch
︡217eefab-b130-4616-a508-d9ff2e8cb623︡{"done":true}
︠343d41bf-fd4a-452e-8e94-8e7731ddc837s︠
fK = f.change_ring(K)
fac1 = list(fK.factor())[2][0] # fac1 modp = x^2
fac2 = list(fK.factor())[0][0] # fac2 modp = (x - 2)^2
fac3 = list(fK.factor())[1][0] # fac3 modp = (x - 3)^2
def coeff(f):
    B = f(0)
    A = f(1) - B - 1
    C = B - A^2/4
    return A, B, C
A1, B1, C1 = coeff(fac1)
A2, B2, C2 = coeff(fac2)
A3, B3, C3 = coeff(fac3)
︡09910b06-42ba-460c-a9f7-e8c37ee72004︡{"done":true}
︠e0200c30-9564-465c-8fe5-5bb6fe9d848es︠
N = 32 # truncation level
rang = 5
TT.<t1> = PowerSeriesRing(K, 't1', default_prec = N)
l1 = t1*(1+C1*t1^2)^(-1/2) # be careful with t1: t1 --> 1/(x + A1/2)
TT.<t2> = PowerSeriesRing(K, 't2', default_prec = N)
l2 = t2*(1+C2*t2^2)^(-1/2) # be careful with t2: t2 --> 1/(x + A2/2)
TT.<t3> = PowerSeriesRing(K, 't3', default_prec = N)
l3 = t3*(1+C3*t3^2)^(-1/2) # be careful with t3: t3 --> 1/(x + A3/2)
︡df07c512-f954-48ee-8ddc-6ede1c88bda2︡{"done":true}
︠883acb1b-c8b7-4778-8aae-37f16df6bd40s︠
# 0) The components v_plus and v_minus
TT.<t1,t2,t3> = PowerSeriesRing(K, 't1,t2,t3', default_prec = 3*N)
l = TT(l1)*TT(l2)*TT(l3)
︡a6a8101c-fdc3-4af4-837a-7bdfe59890b2︡{"done":true}
︠7c33fc25-1334-4630-af78-c7bf11ac85fds︠
# first pole reduction : x-coordinate is -A1/2
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
w = [l.truncate()(1/t,1/(t-A1/2+A2/2),1/(t-A1/2+A3/2))*(t-A1/2)^i/2 for i in range(rang)] # original forms around the first pole
d1 = [w[i].residue() for i in range(rang)] # coefficients of dx/(x+A1/2)
w = [w[i] - d1[i]*(1/t) for i in range(rang)]
ED = [0,0] + [-i/(t^(i+1)) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u> = QQ[]
F1 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F1[i] += FF[i][j]/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A1/2
︡150efe13-e7b2-4c75-8193-c1ff331b8ad6︡{"done":true}
︠85bbac72-9caa-49cf-a907-2464d83eb821s︠
# test
[(l.truncate()(1/t,1/(t-A1/2+A2/2),1/(t-A1/2+A3/2))*(t-A1/2)^i/2 - F1[i](t).derivative() - d1[i]*(1/t)).valuation() for i in range(rang)]
︡d0cd44a5-a557-40fb-a3fc-cd1b9bf51e48︡{"stdout":"[0, 0, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠9c619c10-8713-4759-922d-3b9b6b36819ds︠
# second pole reduction : x-coordinate is -A2/2
w = [l.truncate()(1/(t-A2/2+A1/2),1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i/2 - F1[i](t-A2/2+A1/2).derivative() - d1[i]*(1/(t-A2/2+A1/2)) for i in range(rang)]
d2 = [w[i].residue() for i in range(rang)] # coefficients of dx/(x+A2/2)
w = [w[i] - d2[i]*(1/t) for i in range(rang)]
ED = [0,0] + [-i/(t^(i+1)) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u> = QQ[]
F2 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F2[i] += FF[i][j]/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A2/2
︡82c868e2-9abf-4457-ab0c-94afee1945a0︡{"done":true}
︠34b962d4-968c-468c-b4df-a58881d4bc24s︠
# test
[(l.truncate()(1/(t-A2/2+A1/2),1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i/2 - F1[i](t-A2/2+A1/2).derivative() - d1[i]*(1/(t-A2/2+A1/2)) - F2[i](t).derivative() - d2[i]*(1/t)).valuation()  for i in range(rang)]
︡649512ab-d6f2-42dc-b7b6-38ddc90ed29a︡{"stdout":"[0, 0, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠30cac452-3c7b-4f1b-bd12-eae028e7424es︠
# third pole reduction : x-coordinate is -A3/2
w = [l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i/2 - F1[i](t-A3/2+A1/2).derivative() - d1[i]*(1/(t-A3/2+A1/2)) - F2[i](t-A3/2+A2/2).derivative() - d2[i]*(1/(t-A3/2+A2/2)) for i in range(rang)]
d3 = [w[i].residue() for i in range(rang)] # coefficients of dx/(x+A3/2)
w = [w[i] - d3[i]*(1/t) for i in range(rang)]
ED = [0,0] + [-i/(t^(i+1)) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u> = QQ[]
F3 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F3[i] += FF[i][j]/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A3/2
︡4023f018-2eea-43ed-b0b2-76470ef75ec5︡{"done":true}
︠d0a0c522-3e27-4b35-aa67-7278710d83b5s︠
# test
[(l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i/2 - F1[i](t-A3/2+A1/2).derivative() - d1[i]*(1/(t-A3/2+A1/2)) - F2[i](t-A3/2+A2/2).derivative() - d2[i]*(1/(t-A3/2+A2/2)) - F3[i](t).derivative() - d3[i]*(1/t)).valuation() for i in range(rang)]
︡64e534df-e568-4e05-b304-45b9223aec09︡{"stdout":"[+Infinity, +Infinity, +Infinity, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠3f9ac3dd-f43d-4086-8f82-d48539fcc1cfs︠
dinf30 = (l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3/2 - F1[3](t-A3/2+A1/2).derivative() - d1[3]*(1/(t-A3/2+A1/2)) - F2[3](t-A3/2+A2/2).derivative() - d2[3]*(1/(t-A3/2+A2/2)) - F3[3](t).derivative() - d3[3]*(1/t)).list()[0]
dinf41 = (l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4/2 - F1[4](t-A3/2+A1/2).derivative() - d1[4]*(1/(t-A3/2+A1/2)) - F2[4](t-A3/2+A2/2).derivative() - d2[4]*(1/(t-A3/2+A2/2)) - F3[4](t).derivative() - d3[4]*(1/t)).list()[1]
dinf40 = (l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4/2 - F1[4](t-A3/2+A1/2).derivative() - d1[4]*(1/(t-A3/2+A1/2)) - F2[4](t-A3/2+A2/2).derivative() - d2[4]*(1/(t-A3/2+A2/2)) - F3[4](t).derivative() - d3[4]*(1/t)).list()[0] + dinf41*A3/2
︡68d7c81d-5920-473d-b58b-dd8eaa5785a3︡{"done":true}
︠001d4a0b-c77b-40a5-902c-db8b7ebd739cs︠
# test
[l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i/2 == F1[i](t-A3/2+A1/2).derivative() + d1[i]*(1/(t-A3/2+A1/2)) + F2[i](t-A3/2+A2/2).derivative() + d2[i]*(1/(t-A3/2+A2/2)) + F3[i](t).derivative() + d3[i]*(1/t) for i in range(3)]
l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3/2 == F1[3](t-A3/2+A1/2).derivative() + d1[3]*(1/(t-A3/2+A1/2)) + F2[3](t-A3/2+A2/2).derivative() + d2[3]*(1/(t-A3/2+A2/2)) + F3[3](t).derivative() + d3[3]*(1/t) + dinf30
l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4/2 == F1[4](t-A3/2+A1/2).derivative() + d1[4]*(1/(t-A3/2+A1/2)) + F2[4](t-A3/2+A2/2).derivative() + d2[4]*(1/(t-A3/2+A2/2)) + F3[4](t).derivative() + d3[4]*(1/t) + dinf40 + dinf41*(t-A3/2)
︡5c77e0d7-0104-45aa-92e3-c8a115a1422c︡{"stdout":"[True, True, True]"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"done":true}
︠d9bd3a94-4231-41b7-92fb-735893b23e99s︠
# RESULT:
# for i = 0,1,2, omega_i = dF1[i] + d1[i]*dx/(x+A1/2) + dF2[i] + d2[i]*dx/(x+A2/2) + dF3[i] + d3[i]*dx/(x+A3/2)
# omega_3 = dF1[3] + d1[3]*dx/(x+A1/2) + dF2[3] + d2[3]*dx/(x+A2/2) + dF3[3] + d3[3]*dx/(x+A3/2) + dinf30*dx
# omega_4 = dF1[4] + d1[4]*dx/(x+A1/2) + dF2[4] + d2[4]*dx/(x+A2/2) + dF3[4] + d3[4]*dx/(x+A3/2) + dinf40*dx + dinf41*x*dx
︡410d9eee-30a8-4dc5-801b-04b3fd1d6fac︡{"done":true}
︠b605bf06-8c1e-4233-a8ec-7fe66f0b18eas︠
# 1) The component v_1
TT.<t2,t3> = PowerSeriesRing(K, 't2,t3', default_prec = 3*N)
l = TT(l2)*TT(l3)
g = fac1
gder = g.derivative()
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
︡6ba028f7-ed90-49ae-b33c-d5e24428497b︡{"done":true}
︠63811049-f5ca-47f3-b980-514f2d10c5d8s︠
# first pole reduction : x-coordinate is -A2/2
x = t - A2/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i*x.derivative()/(2*y) for i in range(rang)] # original forms around the first pole
d12 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A2/2)*dx/2y
w = [w[i] - d12[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F12 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F12[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A2/2
︡34ff583f-2d7b-4fc0-a70f-13eec7e4af4f︡{"done":true}
︠d4dafe3e-c25e-4ff8-80f6-9612dfefc33fs︠
# test
[(l.truncate()(1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i*x.derivative()/(2*y) - F12[i](t,y).derivative() - d12[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
︡be1235ea-52eb-4e69-90e9-c214eee363de︡{"stdout":"[0, 0, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠7fc644b0-199c-4b97-a9ba-8365b3150017s︠
# second pole reduction : x-coordinate is -A3/2
x = t - A3/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F12[i](t-A3/2+A2/2,y).derivative() - d12[i]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) for i in range(rang)]
d13 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A3/2)*dx/2y
w = [w[i] - d13[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F13 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F13[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A3/2
︡f460edd2-09d6-4eaa-a2ed-639dc177d818︡{"done":true}
︠10aad194-cfb5-451e-be67-348caafa65e1s︠
# test
[(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F12[i](t-A3/2+A2/2,y).derivative() - d12[i]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[i](t,y).derivative() - d13[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
︡2c9da069-82aa-4c17-94c8-2605388209ff︡{"stdout":"[+Infinity, +Infinity, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠4b5630e5-48ce-4e08-b0e0-a492ad81fcdbs︠
dinf120 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) - F12[2](t-A3/2+A2/2,y).derivative() - d12[2]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[2](t,y).derivative() - d13[2]*(1/t*x.derivative()/(2*y)))).list()[0]
dinf131 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F12[3](t-A3/2+A2/2,y).derivative() - d12[3]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[3](t,y).derivative() - d13[3]*(1/t*x.derivative()/(2*y)))).list()[1]
dinf130 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F12[3](t-A3/2+A2/2,y).derivative() - d12[3]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[3](t,y).derivative() - d13[3]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf131*A3/2
dinf142 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F12[4](t-A3/2+A2/2,y).derivative() - d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[4](t,y).derivative() - d13[4]*(1/t*x.derivative()/(2*y)))).list()[2]
dinf141 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F12[4](t-A3/2+A2/2,y).derivative() - d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[4](t,y).derivative() - d13[4]*(1/t*x.derivative()/(2*y)))).list()[1] + dinf142*A3
dinf140 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F12[4](t-A3/2+A2/2,y).derivative() - d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[4](t,y).derivative() - d13[4]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf141*A3/2 - dinf142*A3^2/4
︡38fc5e37-a426-442b-b025-f53cdb2a04d2︡{"done":true}
︠e8e991ac-b25d-4d22-83e5-3f9ae42cd866s︠
# test
[l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) == F12[i](t-A3/2+A2/2,y).derivative() + d12[i]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[i](t,y).derivative() + d13[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) == F12[2](t-A3/2+A2/2,y).derivative() + d12[2]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[2](t,y).derivative() + d13[2]*(1/t*x.derivative()/(2*y)) + dinf120*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) == F12[3](t-A3/2+A2/2,y).derivative() + d12[3]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[3](t,y).derivative() + d13[3]*(1/t*x.derivative()/(2*y)) + dinf130*x.derivative()/(2*y) + dinf131*x*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) == F12[4](t-A3/2+A2/2,y).derivative() + d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[4](t,y).derivative() + d13[4]*(1/t*x.derivative()/(2*y)) + dinf140*x.derivative()/(2*y) + dinf141*x*x.derivative()/(2*y) + dinf142*x^2*x.derivative()/(2*y)
︡e93ad613-610a-461c-97b9-faa68ff3a25f︡{"stdout":"[True, True]"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"done":true}
︠78b6d945-b337-4dea-8d0b-aeeb42495d50s︠
# RESULT:
# For i = 0,1, omega_i = dF12[i] + d12[i]*1/(x+A2/2)*dx/2y + dF13[i] + d13[i]*1/(x+A3/2)*dx/2y
# omega_2 = dF12[2] + d12[2]*1/(x+A2/2)*dx/2y + dF13[2] + d13[2]*1/(x+A3/2)*dx/2y + dinf120*dx/2y
# omega_3 = dF12[3] + d12[3]*1/(x+A2/2)*dx/2y + dF13[3] + d13[3]*1/(x+A3/2)*dx/2y + dinf130*dx/2y + dinf131*x*dx/2y
# omega_4 = dF12[4] + d12[4]*1/(x+A2/2)*dx/2y + dF13[4] + d13[4]*1/(x+A3/2)*dx/2y + dinf140*dx/2y + dinf141*x*dx/2y + dinf142*x^2*dx/2y
︡571ba0aa-697f-4669-b6f5-f06e7e1bb63a︡{"done":true}
︠6286be27-a15e-4176-8e00-e6bd83ca3d7bs︠
# 2) The component v_2
TT.<t1,t3> = PowerSeriesRing(K, 't1,t3', default_prec = 3*N)
l = TT(l1)*TT(l3)
g = fac2
gder = g.derivative()
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
︡614f5a57-27a1-418b-aa3c-fb5ee4965c75︡{"done":true}
︠be4767ed-c90a-48e5-9ead-7007ba100a7ds︠
# first pole reduction : x-coordinate is -A1/2
x = t - A1/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/t,1/(t-A1/2+A3/2))*(t-A1/2)^i*x.derivative()/(2*y) for i in range(rang)] # original forms around the first pole
d21 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A1/2)*dx/2y
w = [w[i] - d21[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F21 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F21[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A1/2
︡cf27cc6c-3785-4463-b433-5676d85d61b6︡{"done":true}
︠ba2894c5-c2fc-4a76-9137-76916e595597s︠
# test
[(l.truncate()(1/t,1/(t-A1/2+A3/2))*(t-A1/2)^i*x.derivative()/(2*y) - F21[i](t,y).derivative() - d21[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
︡de104905-42dd-4dc3-8ec5-7428f82d83dc︡{"stdout":"[0, 0, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠0c19265f-66e6-4d98-a4c7-80fbea081178s︠
# second pole reduction : x-coordinate is -A3/2
x = t - A3/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F21[i](t-A3/2+A1/2,y).derivative() - d21[i]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) for i in range(rang)]
d23 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A3/2)*dx/2y
w = [w[i] - d23[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F23 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F23[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A3/2
︡8c1268a3-68eb-4c88-800a-aa88ea18dc66︡{"done":true}
︠5f8b66ee-d12e-4636-b558-4521a4f24341s︠
# test
[(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F21[i](t-A3/2+A1/2,y).derivative() - d21[i]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[i](t,y).derivative() - d23[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
︡2c1de1e7-95f9-41ec-bbbe-e18fff2ec753︡{"stdout":"[+Infinity, +Infinity, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠986057bd-f558-40ee-824a-96f3eaba57cas︠
dinf220 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) - F21[2](t-A3/2+A1/2,y).derivative() - d21[2]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[2](t,y).derivative() - d23[2]*(1/t*x.derivative()/(2*y)))).list()[0]
dinf231 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F21[3](t-A3/2+A1/2,y).derivative() - d21[3]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[3](t,y).derivative() - d23[3]*(1/t*x.derivative()/(2*y)))).list()[1]
dinf230 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F21[3](t-A3/2+A1/2,y).derivative() - d21[3]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[3](t,y).derivative() - d23[3]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf231*A3/2
dinf242 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F21[4](t-A3/2+A1/2,y).derivative() - d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[4](t,y).derivative() - d23[4]*(1/t*x.derivative()/(2*y)))).list()[2]
dinf241 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F21[4](t-A3/2+A1/2,y).derivative() - d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[4](t,y).derivative() - d23[4]*(1/t*x.derivative()/(2*y)))).list()[1] + dinf242*A3
dinf240 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F21[4](t-A3/2+A1/2,y).derivative() - d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[4](t,y).derivative() - d23[4]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf241*A3/2 - dinf242*A3^2/4
︡ff905aae-827d-4fea-9eb5-39c66fb9a9e1︡{"done":true}
︠c2e8fe3b-d738-4b22-8097-654f814efd99s︠
# test
[l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) == F21[i](t-A3/2+A1/2,y).derivative() + d21[i]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[i](t,y).derivative() + d23[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) == F21[2](t-A3/2+A1/2,y).derivative() + d21[2]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[2](t,y).derivative() + d23[2]*(1/t*x.derivative()/(2*y)) + dinf220*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) == F21[3](t-A3/2+A1/2,y).derivative() + d21[3]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[3](t,y).derivative() + d23[3]*(1/t*x.derivative()/(2*y)) + dinf230*x.derivative()/(2*y) + dinf231*x*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) == F21[4](t-A3/2+A1/2,y).derivative() + d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[4](t,y).derivative() + d23[4]*(1/t*x.derivative()/(2*y)) + dinf240*x.derivative()/(2*y) + dinf241*x*x.derivative()/(2*y) + dinf242*x^2*x.derivative()/(2*y)
︡b48ffcdd-a9f7-4676-979b-fbb6b30cc1ac︡{"stdout":"[True, True]"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"done":true}
︠8135a9a0-2bf7-4f6d-a875-83d0b1bd431cs︠
# RESULT:
# For i = 0,1, omega_i = dF21[i] + d21[i]*1/(x+A1/2)*dx/2y + dF23[i] + d23[i]*1/(x+A3/2)*dx/2y
# omega_2 = dF21[2] + d21[2]*1/(x+A1/2)*dx/2y + dF23[2] + d23[2]*1/(x+A3/2)*dx/2y + dinf220*dx/2y
# omega_3 = dF21[3] + d21[3]*1/(x+A1/2)*dx/2y + dF23[3] + d23[3]*1/(x+A3/2)*dx/2y + dinf230*dx/2y + dinf231*x*dx/2y
# omega_4 = dF21[4] + d21[4]*1/(x+A1/2)*dx/2y + dF23[4] + d23[4]*1/(x+A3/2)*dx/2y + dinf240*dx/2y + dinf241*x*dx/2y + dinf242*x^2*dx/2y
︡04a92423-c478-4473-a351-93eb4d4d9a92︡{"done":true}
︠f203b4d6-db17-4e3a-a807-b6c9c04eb831s︠
# 3) The component v_3
TT.<t1,t2> = PowerSeriesRing(K, 't1,t2', default_prec = 3*N)
l = TT(l1)*TT(l2)
g = fac3
gder = g.derivative()
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
︡51e407e2-072a-4550-9284-4dcd1beee2be︡{"done":true}
︠f86eecfa-554f-4e9c-a5c6-356df25a41fbs︠
# first pole reduction : x-coordinate is -A1/2
x = t - A1/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/t,1/(t-A1/2+A2/2))*(t-A1/2)^i*x.derivative()/(2*y) for i in range(rang)] # original forms around the first pole
d31 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A1/2)*dx/2y
w = [w[i] - d31[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F31 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F31[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A1/2
︡f4d9eee2-a37a-4a18-a906-aa9dfc3bf20f︡{"done":true}
︠c3780d94-37fa-4572-ac53-6db8ac4e0b37s︠
# test
[(l.truncate()(1/t,1/(t-A1/2+A2/2))*(t-A1/2)^i*x.derivative()/(2*y) - F31[i](t,y).derivative() - d31[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
︡9f7d7fd7-d1ea-49f5-a39b-612144d5d938︡{"stdout":"[0, 0, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠87c034ad-bd72-450c-ac90-e5152d65e6d2s︠
# second pole reduction : x-coordinate is -A2/2
x = t - A2/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^i*x.derivative()/(2*y) - F31[i](t-A2/2+A1/2,y).derivative() - d31[i]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) for i in range(rang)]
d32 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A2/2)*dx/2y
w = [w[i] - d32[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F32 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F32[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A2/2
︡bd0906d1-3d93-46ab-babb-7f5e71173f2c︡{"done":true}
︠a584cbf0-d576-419e-bfe2-3e0201e02c68s︠
# test
[(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^i*x.derivative()/(2*y) - F31[i](t-A2/2+A1/2,y).derivative() - d31[i]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[i](t,y).derivative() - d32[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
︡18095cc2-7eba-4be7-b593-b194d3e209ed︡{"stdout":"[+Infinity, +Infinity, 0, 0, 0]"}︡{"stdout":"\n"}︡{"done":true}
︠0a6f3836-d383-4948-a02b-2c3bde7e2113s︠
dinf320 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^2*x.derivative()/(2*y) - F31[2](t-A2/2+A1/2,y).derivative() - d31[2]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[2](t,y).derivative() - d32[2]*(1/t*x.derivative()/(2*y)))).list()[0]
dinf331 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^3*x.derivative()/(2*y) - F31[3](t-A2/2+A1/2,y).derivative() - d31[3]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[3](t,y).derivative() - d32[3]*(1/t*x.derivative()/(2*y)))).list()[1]
dinf330 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^3*x.derivative()/(2*y) - F31[3](t-A2/2+A1/2,y).derivative() - d31[3]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[3](t,y).derivative() - d32[3]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf331*A2/2
dinf342 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) - F31[4](t-A2/2+A1/2,y).derivative() - d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[4](t,y).derivative() - d32[4]*(1/t*x.derivative()/(2*y)))).list()[2]
dinf341 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) - F31[4](t-A2/2+A1/2,y).derivative() - d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[4](t,y).derivative() - d32[4]*(1/t*x.derivative()/(2*y)))).list()[1] + dinf342*A2
dinf340 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) - F31[4](t-A2/2+A1/2,y).derivative() - d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[4](t,y).derivative() - d32[4]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf341*A2/2 - dinf342*A2^2/4
︡180d2c0a-86a7-4cfc-bbf1-73b1fe99eb99︡{"done":true}
︠21b620d1-0348-4de2-b4e8-ec2ac066a656s︠
# test
[l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^i*x.derivative()/(2*y) == F31[i](t-A2/2+A1/2,y).derivative() + d31[i]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[i](t,y).derivative() + d32[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^2*x.derivative()/(2*y) == F31[2](t-A2/2+A1/2,y).derivative() + d31[2]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[2](t,y).derivative() + d32[2]*(1/t*x.derivative()/(2*y)) + dinf320*x.derivative()/(2*y)
l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^3*x.derivative()/(2*y) == F31[3](t-A2/2+A1/2,y).derivative() + d31[3]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[3](t,y).derivative() + d32[3]*(1/t*x.derivative()/(2*y)) + dinf330*x.derivative()/(2*y) + dinf331*x*x.derivative()/(2*y)
l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) == F31[4](t-A2/2+A1/2,y).derivative() + d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[4](t,y).derivative() + d32[4]*(1/t*x.derivative()/(2*y)) + dinf340*x.derivative()/(2*y) + dinf341*x*x.derivative()/(2*y) + dinf342*x^2*x.derivative()/(2*y)
︡0daab9fa-1fcc-4206-9778-d931fe70d62f︡{"stdout":"[True, True]"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"stdout":"True"}︡{"stdout":"\n"}︡{"done":true}
︠b542205a-1750-43f9-889f-ebd9e26c83c9s︠
# RESULT:
# For i = 0,1, omega_i = dF31[i] + d31[i]*1/(x+A1/2)*dx/2y + dF32[i] + d32[i]*1/(x+A2/2)*dx/2y
# omega_2 = dF31[2] + d31[2]*1/(x+A1/2)*dx/2y + dF32[2] + d32[2]*1/(x+A2/2)*dx/2y + dinf320*dx/2y
# omega_3 = dF31[3] + d31[3]*1/(x+A1/2)*dx/2y + dF32[3] + d32[3]*1/(x+A2/2)*dx/2y + dinf330*dx/2y + dinf331*x*dx/2y
# omega_4 = dF31[4] + d31[4]*1/(x+A1/2)*dx/2y + dF32[4] + d32[4]*1/(x+A2/2)*dx/2y + dinf340*dx/2y + dinf341*x*dx/2y + dinf342*x^2*dx/2y
︡10e31f70-6b7e-4216-bac4-cdcbc4677517︡{"done":true}
︠bfbe40c2-b6e1-4274-b718-82c2bf9b79bds︠
# INTEGRATION

def Int_vplus(i,z1,z2): # integral of omega_i on v_plus from z1 to z2
    x1 = L(z1[0])
    x2 = L(z2[0])
    exact_part1 = F1[i](x2 + A1/2) - F1[i](x1 + A1/2)
    exact_part2 = F2[i](x2 + A2/2) - F2[i](x1 + A2/2)
    exact_part3 = F3[i](x2 + A3/2) - F3[i](x1 + A3/2)
    exact_part = exact_part1 + exact_part2 + exact_part3
    third_kind_part1 = d1[i]*(Log(x2 + A1/2) - Log(x1 + A1/2))
    third_kind_part2 = d2[i]*(Log(x2 + A2/2) - Log(x1 + A2/2))
    third_kind_part3 = d3[i]*(Log(x2 + A3/2) - Log(x1 + A3/2))
    third_kind_part = third_kind_part1 + third_kind_part2 + third_kind_part3
    if i == 3:
        inf_part = dinf30*(x2 - x1)
    elif i == 4:
        inf_part = dinf40*(x2 - x1) + dinf41/2*(x2^2 - x1^2)
    else:
        inf_part = 0
    return exact_part + third_kind_part + inf_part

def Int_vminus(i,z1,z2): # integral of omega_i on v_minus from z1 to z2
    return -Int_vplus(i,z1,z2)

def NC(z,j): # points on v_j in the new coordinates
    x, y = L(z[0]), L(z[1])
    Z1, Z2, Z3 = x + A1/2, x + A2/2, x + A3/2
    l1 = Z1*(1 + C1/Z1^2).sqrt()
    l2 = Z2*(1 + C2/Z2^2).sqrt()
    l3 = Z3*(1 + C3/Z3^2).sqrt()
    if j == 1:
        l = l2*l3
    elif j == 2:
        l = l1*l3
    elif j == 3:
        l = l1*l2
    return (x,y/l)

def Int(i,j,z1,z2): # integral of omega_i on v_j from z1 to z2
    RR.<T> = K[]
    x1, y1 = NC(z1,j)
    x2, y2 = NC(z2,j)
    if j == 1:
        exact_part2 = F12[i](x2+A2/2,y2) - F12[i](x1+A2/2,y1)
        exact_part3 = F13[i](x2+A3/2,y2) - F13[i](x1+A3/2,y1)
        exact_part = exact_part2 + exact_part3
        pol12 = T^2 + (A2-A1)*T - C1
        pol13 = T^2 + (A3-A1)*T - C1
        r12, s12 = pol12.roots(multiplicities=False)
        r13, s13 = pol13.roots(multiplicities=False)
        T1 = x1 + y1 + A1/2
        T2 = x2 + y2 + A1/2
        third_kind_part_2 = d12[i]/(r12-s12)*(Log((T2-r12)/(T2-s12)) - Log((T1-r12)/(T1-s12)))
        third_kind_part_3 = d13[i]/(r13-s13)*(Log((T2-r13)/(T2-s13)) - Log((T1-r13)/(T1-s13)))
        third_kind_part = third_kind_part_2 + third_kind_part_3
        if i == 2:
            inf_part0 = dinf120/2*(Log(T2) - Log(T1))
            inf_part1 = 0
            inf_part2 = 0
        elif i == 3:
            inf_part0 = dinf130/2*(Log(T2) - Log(T1))
            inf_part1 = dinf131/4*((T2-A1*Log(T2)+C1/T2) - (T1-A1*Log(T1)+C1/T1))
            inf_part2 = 0
        elif i == 4:
            inf_part0 = dinf140/2*(Log(T2) - Log(T1))
            inf_part1 = dinf141/4*((T2-A1*Log(T2)+C1/T2) - (T1-A1*Log(T1)+C1/T1))
            inf_part2 = dinf142/8*((T2^2/2-2*A1*T2+(A1^2-2*C1)*Log(T2)-2*A1*C1/T2-C1^2/(2*T2^2)) - (T1^2/2-2*A1*T1+(A1^2-2*C1)*Log(T1)-2*A1*C1/T1-C1^2/(2*T1^2)))
        else:
            inf_part0 = 0
            inf_part1 = 0
            inf_part2 = 0
        inf_part = inf_part0 + inf_part1 + inf_part2
    elif j == 2:
        exact_part1 = F21[i](x2+A1/2,y2) - F21[i](x1+A1/2,y1)
        exact_part3 = F23[i](x2+A3/2,y2) - F23[i](x1+A3/2,y1)
        exact_part = exact_part1 + exact_part3
        pol21 = T^2 + (A1-A2)*T - C2
        pol23 = T^2 + (A3-A2)*T - C2
        r21, s21 = pol21.roots(multiplicities=False)
        r23, s23 = pol23.roots(multiplicities=False)
        T1 = x1 + y1 + A2/2
        T2 = x2 + y2 + A2/2
        third_kind_part_1 = d21[i]/(r21-s21)*(Log((T2-r21)/(T2-s21)) - Log((T1-r21)/(T1-s21)))
        third_kind_part_3 = d23[i]/(r23-s23)*(Log((T2-r23)/(T2-s23)) - Log((T1-r23)/(T1-s23)))
        third_kind_part = third_kind_part_1 + third_kind_part_3
        if i == 2:
            inf_part0 = dinf220/2*(Log(T2) - Log(T1))
            inf_part1 = 0
            inf_part2 = 0
        elif i == 3:
            inf_part0 = dinf230/2*(Log(T2) - Log(T1))
            inf_part1 = dinf231/4*((T2-A2*Log(T2)+C2/T2) - (T1-A2*Log(T1)+C2/T1))
            inf_part2 = 0
        elif i == 4:
            inf_part0 = dinf240/2*(Log(T2) - Log(T1))
            inf_part1 = dinf241/4*((T2-A2*Log(T2)+C2/T2) - (T1-A2*Log(T1)+C2/T1))
            inf_part2 = dinf242/8*((T2^2/2-2*A2*T2+(A2^2-2*C2)*Log(T2)-2*A2*C2/T2-C2^2/(2*T2^2)) - (T1^2/2-2*A2*T1+(A2^2-2*C2)*Log(T1)-2*A2*C2/T1-C2^2/(2*T1^2)))
        else:
            inf_part0 = 0
            inf_part1 = 0
            inf_part2 = 0
        inf_part = inf_part0 + inf_part1 + inf_part2
    elif j == 3:
        exact_part1 = F31[i](x2+A1/2,y2) - F31[i](x1+A1/2,y1)
        exact_part2 = F32[i](x2+A2/2,y2) - F32[i](x1+A2/2,y1)
        exact_part = exact_part1 + exact_part2
        pol31 = T^2 + (A1-A3)*T - C3
        pol32 = T^2 + (A2-A3)*T - C3
        r31, s31 = pol31.roots(multiplicities=False)
        r32, s32 = pol32.roots(multiplicities=False)
        T1 = x1 + y1 + A3/2
        T2 = x2 + y2 + A3/2
        third_kind_part_1 = d31[i]/(r31-s31)*(Log((T2-r31)/(T2-s31)) - Log((T1-r31)/(T1-s31)))
        third_kind_part_2 = d32[i]/(r32-s32)*(Log((T2-r32)/(T2-s32)) - Log((T1-r32)/(T1-s32)))
        third_kind_part = third_kind_part_1 + third_kind_part_2
        if i == 2:
            inf_part0 = dinf320/2*(Log(T2) - Log(T1))
            inf_part1 = 0
            inf_part2 = 0
        elif i == 3:
            inf_part0 = dinf330/2*(Log(T2) - Log(T1))
            inf_part1 = dinf331/4*((T2-A3*Log(T2)+C3/T2) - (T1-A3*Log(T1)+C3/T1))
            inf_part2 = 0
        elif i == 4:
            inf_part0 = dinf340/2*(Log(T2) - Log(T1))
            inf_part1 = dinf341/4*((T2-A3*Log(T2)+C3/T2) - (T1-A3*Log(T1)+C3/T1))
            inf_part2 = dinf342/8*((T2^2/2-2*A3*T2+(A3^2-2*C3)*Log(T2)-2*A3*C3/T2-C3^2/(2*T2^2)) - (T1^2/2-2*A3*T1+(A3^2-2*C3)*Log(T1)-2*A3*C3/T1-C3^2/(2*T1^2)))
        else:
            inf_part0 = 0
            inf_part1 = 0
            inf_part2 = 0
        inf_part = inf_part0 + inf_part1 + inf_part2
    return exact_part + third_kind_part + inf_part
︡b0e5ddc6-fdff-4ad9-a6f3-792c3cdcdebe︡{"done":true}
︠3d0c1d56-6a43-46f0-854f-74c9557bfb52s︠
# period integrals
per1 = [Int(i,1,P1,P2) + Int_vplus(i,P2,P3) + Int(i,2,P3,P4) + Int_vminus(i,P4,P1) for i in range(rang)]
per2 = [Int(i,2,P3,P4) + Int_vminus(i,P4,P5) + Int(i,3,P5,P6) + Int_vplus(i,P6,P3) for i in range(rang)]
︡25839b83-1e59-4147-8fab-f64474164cea︡{"done":true}
︠eca7d0d2-0afd-447c-bf7d-9032a07378c9s︠
path = [Int_vminus(i,S,P1) + Int(i,1,P1,P2) + Int_vplus(i,P2,R) for i in range(rang)] # image of path under tau = e1e2
for i in range(5):
    path[i] - 2/3*per1[i] + 1/3*per2[i]
︡49464027-cf3b-4c33-8ecc-b672603774a5︡{"stdout":"2*a^34 + a^35 + a^38 + a^39 + 4*a^41 + 4*a^43 + O(a^44)\n2*a^34 + a^35 + a^37 + a^38 + 3*a^39 + 3*a^40 + a^42 + a^43 + O(a^44)\na^34 + 3*a^35 + 2*a^36 + 2*a^38 + a^40 + 3*a^41 + 4*a^42 + a^43 + O(a^44)\n1 + 3*a^4 + a^8 + 3*a^12 + a^16 + 3*a^20 + a^24 + 3*a^28 + a^32 + 3*a^34 + 4*a^35 + 3*a^36 + a^37 + 4*a^38 + a^39 + 4*a^40 + 2*a^41 + a^42 + 4*a^43 + O(a^44)\na^34 + 2*a^35 + 3*a^36 + 2*a^38 + 4*a^39 + a^40 + a^41 + 2*a^42 + a^43 + O(a^44)\n"}︡{"done":true}
︠93a2e63c-a2cd-4561-82bf-72c558c57ac3︠
path = [Int_vminus(i,S,P4) + Int(i,2,P4,P3) + Int_vplus(i,P3,R) for i in range(rang)] # image of path under tau = (-e4)(-e3)
for i in range(5):
    path[i] + 1/3*per1[i] + 1/3*per2[i]
︡89cfa44f-4224-45c2-90bb-2c62b8b2ba52︡{"stdout":"2*a^34 + a^35 + a^38 + a^39 + 4*a^41 + 4*a^43 + O(a^44)\n2*a^34 + a^35 + a^37 + a^38 + 3*a^39 + 3*a^40 + a^42 + a^43 + O(a^44)\na^34 + 3*a^35 + 2*a^36 + 2*a^38 + a^40 + 3*a^41 + 4*a^42 + a^43 + O(a^44)\n1 + 3*a^4 + a^8 + 3*a^12 + a^16 + 3*a^20 + a^24 + 3*a^28 + a^32 + 3*a^34 + 4*a^35 + 3*a^36 + a^37 + 4*a^38 + a^39 + 4*a^40 + 2*a^41 + a^42 + 4*a^43 + O(a^44)\na^34 + 2*a^35 + 3*a^36 + 2*a^38 + 4*a^39 + a^40 + a^41 + 2*a^42 + a^43 + O(a^44)\n"}︡{"done":true}
︠d17f6710-6db1-4fb2-b46f-b82e8d68b248︠
path = [Int_vminus(i,S,P5) + Int(i,3,P5,P6) + Int_vplus(i,P6,R) for i in range(rang)] # image of path under tau = e5e6
for i in range(5):
    path[i] + 1/3*per1[i] - 2/3*per2[i]
︡777ed465-8afc-4811-9d41-dd8fe16334df︡{"stdout":"2*a^34 + a^35 + a^38 + a^39 + 4*a^41 + 4*a^43 + O(a^44)\n2*a^34 + a^35 + a^37 + a^38 + 3*a^39 + 3*a^40 + a^42 + a^43 + O(a^44)\na^34 + 3*a^35 + 2*a^36 + 2*a^38 + a^40 + 3*a^41 + 4*a^42 + a^43 + O(a^44)\n1 + 3*a^4 + a^8 + 3*a^12 + a^16 + 3*a^20 + a^24 + 3*a^28 + a^32 + 3*a^34 + 4*a^35 + 3*a^36 + a^37 + 4*a^38 + a^39 + 4*a^40 + 2*a^41 + a^42 + 4*a^43 + O(a^44)\na^34 + 2*a^35 + 3*a^36 + 2*a^38 + 4*a^39 + a^40 + a^41 + 2*a^42 + a^43 + O(a^44)\n"}︡{"done":true}
︠17b0c7f4-64bf-4601-a914-06d521ef1516︠









