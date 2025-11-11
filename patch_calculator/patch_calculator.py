import scipy.integrate as integrate
import scipy.special as special
from numpy import sin, cos, pi, sqrt, arccos

# Constantes
c = 299792458e3  # Vitesse de la lumière dans le vide (mm/s)
Er = 2.2  # Permittivité relative du substrat

f = 2.235e9   # Fréquence de fonctionnement (Hz)

k0 = 2 * pi * f / c  # Nombre d'onde dans le vide (rad/mm)
mu0 = 4 * pi * 1e-7 * 1e6  # Perméabilité du vide (H/mm)
epsilon0 = 8.854e-12 * 1e3  # Permittivité du vide (F/mm)
h = 0.508  # Épaisseur du substrat (mm)

# Pour une antenne patch rectangulaire standard
def largeur_patch(f, Er, c):
    return c/(2*f)*sqrt(2/(Er+1)) #en mm

W = largeur_patch(f, Er, c)
Eff = (Er + 1)/2 + (Er - 1)/(2*sqrt(1 + 12*(h/W))) # Permittivité effective

def deltaL(Eff, W, h):
    return 0.412*h*(Eff+0.3)*(W/h+0.264)/((Eff-0.258)*(W/h+0.8)) #en mm

def longueur_patch(f, Eff, c):
    return ((c/(2.235e9*2*sqrt(2))) - 2*deltaL(Eff, W, h)) #en mm

def Leff(L, deltaL):
    return L + 2*deltaL #en mm






#Pour une antenne patch rectangulaire avec insertion de fentes
X = k0*W

def Si(X):
    return special.sici(k0*W)[0] #sinus integral

print(Si(X))
I1 = -2*cos(X)+X*Si(X)+sin(X)/X

print("I1 : ",I1)

G1 = I1/(120*pi**2) # en S

print("G1 : ", G1)


Rin = 50 #en Ohm : impédance d'entrée

def y0 (Rin, G1, G12, L, W, k0):
    return ((2/pi)*arccos(Rin*2*(G1 + G12(k0, W, L))) )

def G12 (k0, W, L):
    integrale_inter = integrate.quad(lambda theta: (sin(k0*W/2*cos(theta))/cos(theta))**2*special.jv(0, k0*L*sin(theta))*sin(theta)**3, 0, pi)
    return (1/(120*pi**2))*integrale_inter[0] #on ne tient pas compte de l'erreur de l'intégration (elle est négligeable devant la valeur de l'intégrale)

print(G12(46.842e-3, 53.022, 44.161))

# Dimensions d'une antenne patch "basique" pour f = 2235 MHz = 2.235 GHz
print("Largeur du patch :", largeur_patch(f, Er, c))
print("Delta L:", deltaL(Eff, W, h))
print("Longueur du patch :", longueur_patch(f, Eff, c))

# Calcul de la profondeur de l'insert
print("Insert:", y0(Rin, G1, G12(k0, largeur_patch(f, Er, c), longueur_patch(f, Eff, c))))
