def read_params(filename):
   """Read Parameters from file and return Params (P) object"""
   P = Params()
   with open(filename) as f:
       for line in f:
           tmp = [x.strip() for x in line.split("=")]
           sym = tmp[0]
           val = [x.strip() for x in tmp[1].split()][0]
           setattr(P, sym, eval(val))
   return P

class Params(object):
   pass

if __name__ == "__main__":
    P = read_params("params.txt")
    for k in P.__dict__:
       print "{} : {}".format(k, P.__dict__[k])
