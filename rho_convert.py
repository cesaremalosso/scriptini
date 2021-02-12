#!/usr/bin/python
from scipy.constants import N_A

def len(hmass,omass,nmol,rho):
    V = (omass*nmol+hmass*2*nmol)/N_A/rho
    print(str(V**(1/3)*1e9) + ' angstrom')




def main():
    import argparse
    parser = argparse.ArgumentParser(description='Evaluate box simulation length starting from the system density')
    parser.add_argument('-n','--nmols',
                        type=int,
                        required=True,
                        help='number of water molecules')
    parser.add_argument('-t', '--target',
                        type=float,
                        required=True,
                        help='target density')
    parser.add_argument('-m','--masses',
                        type = float,
                        nargs='+',
                        required=True,
                        help='mass of hydrogen and oxygen atoms')
    args = parser.parse_args()
    len(args.masses[0],args.masses[1],args.nmols,args.target)
if __name__ == "__main__":
  main()

