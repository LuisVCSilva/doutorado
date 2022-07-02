from numpy import *
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
from tqdm import tqdm
import json
import sys


global epsilon
epsilon = float("nan")

global hbar2m
hbar2m = float("nan")

global sigma
sigma = float("nan")

def pts_iniciais (r, u):
	u[0] = exp(-sqrt(epsilon/hbar2m*(sigma**12)/25.0)*(r[0]**-5))
	u[1] = exp(-sqrt(epsilon/hbar2m*(sigma**12)/25.0)*(r[1]**-5)) 
	return u

def main(args):
   params = open(args[0],"r").read()
   params = json.loads(params)
   
   min_x = float(params["min_x"])
   max_x = float(params["max_x"])
   malha = int(params["malha"])
   min_energia = float(params["min_energia"])
   max_energia = float(params["max_energia"])
   delta_e = float(params["delta_e"])
   max_l = int(params["max_l"])
   output_file = params["output_file"]
   
   global epsilon
   epsilon = float(params["epsilon"])
   
   global hbar2m   
   hbar2m = float(params["hbar2m"])
   
   global sigma   
   sigma = float(params["sigma"])

   if output_file=="":
      print("")
   else:
      sys.stdout = open(output_file, "w")


   e = min_energia + delta_e
   q = sqrt(e/hbar2m)
   lambda2 = pi/q/2.0
   rmin = min_x
   max_r = max_x + lambda2

   vpot	= zeros(malha+1)
   secao = zeros(max_l+1)

   r = zeros(malha+1)
   u = zeros(malha+1)
   p = zeros(malha+1)
   f = zeros(malha+1)

   ne = int(((max_energia-min_energia)/delta_e) + 1.5)
   #print(ne)
   k = 0
   for n in (range(ne+1)[1:]):  
      e = min_energia + (n-1)*delta_e
      q = sqrt(e/hbar2m)
      lambda2 = pi/q/2.0

      for l in range(max_l+1):
         rmin	= min_x
         secao[l] = 0.0

         max_r = max_x + lambda2
         dx = (max_r-rmin)/malha
         ddx12 = dx*dx/12.0
         r2 = max_r
         i2 = malha
         r1 = max_r - dx* int(lambda2/dx + 0.5)
         i1 = malha - int( lambda2/dx + 0.5 )

         r = [rmin + (i+0.0) * dx for i in range(malha+1)]
         vpot = [epsilon * (-2.0*pow(sigma/r_i,6) + (sigma/r_i**12)) for r_i in r]
         
         f =  [1.0 - ddx12/hbar2m * ( hbar2m*l*(l+1)/(r_i**2) + vpot_i - e ) for (r_i,vpot_i) in zip(r,vpot)]
         
         u = zeros(malha+1)
         pts_iniciais (r, u)
         for i in range(malha)[1:]:
            u[i+1] = ((12.0-10.0*f[i])*u[i]-f[i-1]*u[i-1])/f[i+1]




         normalizacao = 0.0 
         normalizacao = sum([(u_i**2)*dx for u_i in u])
         u = [u_i/sqrt(normalizacao) for u_i in u]



         kappa = r1*u[i2]/(r2*u[i1])
         tandelta = (kappa*spherical_jn(l,q*r1)-spherical_jn(l,q*r2)) / (kappa*spherical_yn(l,q*r1)-spherical_yn(l,q*r2))
         delta = arctan(tandelta)
         secao[l] = secao[l] + 4*pi/(q*q) * (2*l+1)*(sin(delta)**2)
         p = [sin(q*r_i-l*pi/2.0 + delta) for r_i in r]
         p = [p_i/(p[i2]/u[i2]) for p_i in p]
         p = [sin(q*r_i-l*pi/2.0 + delta) for r_i in r]
         p = [p_i/(p[i2]/u[i2]) for p_i in p]

         #plt.plot(u)
         #plt.savefig(str(k)+".png")
         #plt.clf()
         #k = k+1
			 

      print("{} {}".format(e, sum(secao)),end="")
      for l in range(max_l+1):
         print(" {}".format(secao[l]),end=" ")
      print("")


if __name__ == "__main__":
    main(sys.argv[1:])
