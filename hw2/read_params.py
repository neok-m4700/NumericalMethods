def read_params(filename):
   P = Params()
   with open(filename) as f:
       for line in f:
           tmp = [x.strip() for x in line.split("=")]
           sym = tmp[0]
           val = [x.strip() for x in tmp[1].split()][0]
           setattr(P, sym, eval(val))
   #Set number of steps in each direction
   P.Nx = int(P.Lx/P.dx)    
   P.Ny = int(P.Ly/P.dx)    
   P.Nz = int(P.Lz/P.dx)    
   #Set total number of steps
   P.N = P.Nx * P.Ny * P.Nz
   #Set total number of time steps
   P.nSteps = int(P.T/P.dt)
   return P

class Params(object):
   pass

if __name__ == "__main__":
    P = read_params("params.txt")
    for k in P.__dict__:
       print "{} : {}".format(k, P.__dict__[k])
