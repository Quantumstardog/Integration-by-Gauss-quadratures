import numpy as np
from numpy import *


def Legendre(n,x):
    x=np.array(x)
    if (n==0):
        return x*0 + 1.
    elif (n==1):
        return x
    else:
        return ((2.0*n - 1.0)*x*Legendre(n-1,x) - (n-1)*Legendre(n-2,x))/n


#Derivada de los poli. de Legendre

def DLegendre(n,x):
    x=np.array(x)
    if (n==0):
        return x*0
    elif (n==1):
        return x*0 +1.0
    else:
        return (n/(x**2 -1.0))*(x*Legendre(n,x) - Legendre(n-1,x))

#Se define el Polinomio de Lobatto y su derivada mediante recurrencia

def Lobatto(n,x):
    x=np.array(x)
    return   DLegendre(n,x)
      
def DLobatto(n,x):
    x=np.array(x)
    if (n==0):
         return 0
    elif (n==1):
         return 0
    elif (n==2):
         return 3
    elif (n==3):
         return 15*x
    else:
         return ((2*n - 1)*x*DLobatto(n-1,x) - (n+1)*DLobatto(n-2,x))/(n-2)
#Obtenemos las raíces    
def RadauRoots(polyorder, tolerance=1e-20):
    if polyorder <2:
        global err
        err=1
    else:
        roots=[]
        s=-1
        roots.append(s)
 
        #Los polinomios son alternadamente funciones pares e impares. Solo evaluamos la mitad de raices
        for i in range(2,int((polyorder)/2) + 1 ):
            x=np.cos(np.pi*(i - 0.25  )/(polyorder + 0.5))
            error= 10*tolerance
            iters=0
            while (error> tolerance) and (iters < 1000):
                dx= -Lobatto(polyorder-1,x)/DLobatto(polyorder-1,x)
                x= x+dx
                iters= iters+1
                error=abs(dx)
            roots.append(x)
        #Simetria para obtener las raices faltantes
        roots=np.array(roots)
        if polyorder%2==0:
            roots=np.concatenate((-1.0*roots, roots[::-1]))
        else:
            roots=np.concatenate((-1.0*roots, [0.0], roots[::-1]))        
        roots=np.array(roots)
        roots=np.sort(roots)

        error=0 #se encontraron las raices
    return roots,err    



##################################################################
# Weight coefficients


def GaussRadauWeights(polyorder):
     W=[]
     p=2/(polyorder*(polyorder-1))
     W.append(p)
     xis,err=RadauRoots(polyorder)
     if err==0:
          for i in range (1,polyorder-1):
                d=2.0/((polyorder*(polyorder-1))*(Legendre(polyorder-1,xis[i])**2))
                W.append(d)
          g=2/(polyorder*(polyorder-1))
          W.append(g)
          
          err=0
     else:
          err=1
     return W, xis, err


##################################################################
# The integral value 
# func 		: the integrand
# a, b 		: lower and upper limits of the integral
# polyorder 	: order of the Legendre polynomial to be used
#
def GaussRadauQuadrature(func, polyorder, a, b):
	Ws,xs, err= GaussRadauWeights(polyorder)
	if err==0:
		ans=(b-a)*0.5*sum( Ws*func( (b-a)*0.5*xs+ (b+a)*0.5 ) )
	else: 
		# (in case of error)
		err=1
		ans=None
	return ans,err
##################################################################
# The integrand - change as required
def func(x):
	return (4.0/(1+x**2))
##################################################################
# 
 
order=5
Ws=[]
xs=[]
err=0

Ws,xs,err=GaussRadauWeights(order)

if err==0:
	print ("Order    : ", order)
	print ("Roots    : ", xs)
	print ("Weights  : ", Ws)
else:
	print ("Roots/Weights evaluation failed")
 
# Integrating the function
[ans,err]=GaussRadauQuadrature(func , order, 0,1)
if err==0:
	print ("Integral : ", ans)
else:
	print ("Integral evaluation failed")
