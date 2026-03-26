import numpy as np
import matplotlib.pyplot as plt



phi = 9.936667 * np.pi / 180
Gsc = 1353
cita = 9 * np.pi / 180
rho = 0.3
a = 0.23
b = 0.35

t1 = np.array([9,10,11,12,13,14,15])
t2 = t1 + 1

def calcular_mes(n, epsilon, nr, Ib, phi, cita, Gsc, a, b, rho, t1, t2):

    delta = (23.45*np.sin((360*(284+n))/365))*(np.pi/180)
    ws = np.arccos(-np.tan(delta)*np.tan(phi))
    wsgrad = ws*(180/np.pi)
    N = (2/15)*wsgrad

    Rb = (
        np.sin(delta)*np.sin(phi-cita)*(t2-t1)
        + (12/np.pi)*np.cos(delta)*np.cos(phi-cita)
        * (np.sin(np.deg2rad(15*t1)) - np.sin(np.deg2rad(15*t2)))
    ) / (
        np.sin(delta)*np.sin(phi)*(t2-t1)
        + (12/np.pi)*np.cos(delta)*np.cos(phi)
        * (np.sin(np.deg2rad(15*t1)) - np.sin(np.deg2rad(15*t2)))
    )

    Io = Gsc*epsilon*(
        np.sin(delta)*np.sin(phi)*(t2-t1)
        + (12/np.pi)*np.cos(delta)*np.cos(phi)
        * (np.sin(np.deg2rad(15*t1)) - np.sin(np.deg2rad(15*t2)))
    )

    I = Io*(a+b*(nr/N))
    kt = I/Io

    Id = np.where(kt <= 0.35, I*(1-0.249*kt), I*(1.577-1.84*kt))

    It = Ib*Rb + Id*((1+np.cos(cita))/2) + I*((1-np.cos(cita))/2)*rho

    return It

Ite = calcular_mes(16, 1.0343, 6.7, np.array([554,544,518,495,448,403,354]), phi, cita, Gsc, a, b, rho, t1, t2)

Itf = calcular_mes(45, 1.0256, 7.5, np.array([637,638,612,571,511,462,387]), phi, cita, Gsc, a, b, rho, t1, t2)

Itm = calcular_mes(75, 1.0102, 7.9, np.array([654,661,635,581,508,442,353]), phi, cita, Gsc, a, b, rho, t1, t2)

Ita = calcular_mes(106, 0.9927, 7, np.array([543,538,518,447,362,291,230]), phi, cita, Gsc, a, b, rho, t1, t2)

ItM = calcular_mes(136, 0.9771, 4.5, np.array([479,431,398,356,227,168,126]), phi, cita, Gsc, a, b, rho, t1, t2)

Itj = calcular_mes(166, 0.9682, 3.9, np.array([421,378,311,340,212,140,115]), phi, cita, Gsc, a, b, rho, t1, t2)

ItJ = calcular_mes(197, 0.9673, 3.7, np.array([340,322,279,331,203,149,132]), phi, cita, Gsc, a, b, rho, t1, t2)

ItA = calcular_mes(228, 0.975, 3.8, np.array([421,383,342,258,201,161,139]), phi, cita, Gsc, a, b, rho, t1, t2)

Its = calcular_mes(258, 0.9891, 4, np.array([514,432,337,244,173,138,116]), phi, cita, Gsc, a, b, rho, t1, t2)

Ito = calcular_mes(289, 1.007, 4.2, np.array([430,367,291,214,175,140,130]), phi, cita, Gsc, a, b, rho, t1, t2)

Itn = calcular_mes(319, 1.0231, 4.6, np.array([372,331,289,236,205,179,177]), phi, cita, Gsc, a, b, rho, t1, t2)

Itd = calcular_mes(350, 1.0334, 5.3, np.array([482,464,428,379,349,322,277]), phi, cita, Gsc, a, b, rho, t1, t2)
    
    
    
    #MATRIZ It (todos los meses, todas las horas)
   
ItG = np.array([Ite, Itf, Itm, Ita, ItM, Itj, ItJ, ItA, Its, Ito, Itn, Itd])
   
    #graficos
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

c = np.array([16, 45, 75, 106, 136, 166, 197, 228, 258, 289, 319, 350])
a = np.array([9,10,11,12,13,14,15])

A, C = np.meshgrid(a, c)

ax.plot_surface(A, C, ItG, alpha=0.5)

ax.set_xlabel('Hora')
ax.set_ylabel('Día del año')
ax.set_zlabel('Radiación (Wh/m²)')

plt.savefig("results/Radiacion.png", dpi=300)

plt.show()
   
   
   ################################################
 
 
mvidrio = 1.1*1.114*0.004*2800;
cpvidrio = 750;
Tvidrio= 35.5;
Tambiente= 22;
L=1.114;
Tf = (Tvidrio+Tambiente)/2  ;

#Pr y visco tomados al Tf de 30, tabla A15
Pr= 0.7282;
visco= 1.608*(10**(-5));
k = 0.02588;
emisividadvidrio=0.8;
emisividadagua=0.96;
Boltz = 5.6697*(10**(-8));
Tagua = 55;

Pw = 15758;
Pwv = 5628;

beta = 1/(273+ Tf);

#Qs = mvidrio*cpvidrio * deltaTvidrio;

Ray= (9.78* beta*(Tvidrio-Tambiente)*(L**3)*Pr)/(visco**2);

Nu= (0.825+((0.387*(Ray**(1/6)))/((1+(0.492/Pr)**(9/16))**(8/27))))**2;

hsup = (k*Nu)/L;

H = np.array([64, 62.6, 62.2, 65.2, 72.75, 76.25, 70.75, 73, 79.8, 81.2, 75, 67.6])

Trocio=((H/100)**(1/8))*((110+Tambiente)-110)+273; 

Tcielo= (Tambiente+273)*(0.8+((Trocio-273)/(250)))**0.25;

qvidrio=emisividadvidrio*Boltz*(((Tvidrio+273)**4) - (Tcielo**4));

qconv = hsup*(Tvidrio-Tambiente)

#calores p rdidos agua
qra=(Boltz/((1/emisividadagua)+(1/emisividadvidrio)-1))*(((Tagua+273)**4) - ((Tvidrio+273)**4));


hconvnat= 0.884*((Tagua+273)-(Tvidrio+273)+(((Pw-Pwv)/(268.9*(10**3)-Pw))*(Tagua+273+273)))**(1/3); #W/m2K

hT = hconvnat;
#hT es de la cubierta.



he= 0.013*hT;


qcw= hT*(Tagua-Tvidrio);

qevap = he*(Pw-Pwv);

#balance de calor en la parte superior

tau = 0.86;

#Iwat = It/7; #W
Iwat = ItG;
qtrans= Iwat*tau; #W



#balance de calor en la parte inferior

alfaagua = 0.8;
Abandeja= 1.08*0.7; #m2


qreal = qtrans*alfaagua - ((qra/7)+(qcw/7)+(qevap/7));
Qreal= qreal*Abandeja; #W 

Qrreal= Qreal.T;
Qrrealm = np.sum(Qrreal);
Qrrealmp = np.mean(Qrreal);

Acubierta = 1.1*1.114;

#Calor para calentar el vidrio
Qv = (mvidrio*cpvidrio*(Tvidrio-Tambiente))/(3600*7); #W

volumenagua= 1.08*0.7*0.01; #m3

rhoagua= 997;
cpagua= 4180;

masaagua=volumenagua*rhoagua
Qca= (masaagua*cpagua*(Tagua-Tambiente))/(3600*7); #W

entalpiavap= 1859*10^3 ;

Qevapo= (1*masaagua*entalpiavap)/3600;

Qcv = ((qconv/7)*Acubierta) + Qv;

#Q TOTAL EVAPORANDO TODA EL AGUA ES:
QN100 = Qcv + Qca + Qevapo;

ef= (Qreal/QN100)*100;
efe= ef.T;
efecp = np.mean(efe);

md= ((Qreal- Qca)/(masaagua*entalpiavap))*(3600*7);
mdest= md.T;
mdestm = np.sum(mdest)
#graficos

fig = plt.figure()

ax= fig.add_subplot(111, projection ='3d')
c= np.array([16, 45, 75, 106, 136, 166, 197, 228, 258, 289, 319, 350]);
a = np.array([9, 10, 11, 12, 13, 14, 15])

A, C = np.meshgrid (a, c)

ax.plot_surface (A, C, Qreal, alpha=0.5)

ax.set_ylabel ('Dia del año')
ax.set_xlabel ('Hora militar')
ax.set_zlabel ('Calor neto disponible (W)')
plt.savefig("results/CalorNeto.png", dpi=300)

plt.show()
    
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

A, C = np.meshgrid(a, c)

ax.plot_surface(A, C, md, alpha=0.5)

ax.set_xlabel('Hora')
ax.set_ylabel('Día del año')
ax.set_zlabel('Masa de agua destilada (kg)')

plt.savefig("results/Masadeagua.png", dpi=300)

plt.show()

mdest_mes = np.sum(mdest, axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)

d = ['Enero', 'Febrero', 'Marzo', 'Abril', 'Mayo', 'Junio',
     'Julio', 'Agosto', 'Setiembre', 'Octubre', 'Noviembre', 'Diciembre']

ax.bar(d, mdest_mes)

ax.set_xlabel('Mes')
ax.set_ylabel('Masa de agua destilada (kg)')

plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("results/MasadeAguaDestilada.png", dpi=300)
plt.show()