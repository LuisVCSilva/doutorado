import matplotlib.pyplot as plt
import sys
from numpy import *

def main(args):
   arq = open(args[0],"r").read()
   raw = arq.strip().split("\n")
   raw = [list(filter(None, x.strip().split(" "))) for x in raw]
   raw = array(raw).astype(float)
   lineObjects = plt.plot(raw)
   plt.xlabel("E (meV)$")
   plt.ylabel("Seção eficaz (barn)")
   plt.legend(iter(lineObjects), ("Seção de choque $E$", "$\\sigma(E)$", "$\\sigma_l(E)$", "$l=0$", "$l_{max}$"))
   plt.savefig("resultado.png")
   #plt.close()


if __name__ == "__main__":
    main(sys.argv[1:])  
