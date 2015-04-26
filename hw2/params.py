class Params(object):
   def __init__(self,filename):
      """Read Parameters from file and return Params (P) object"""
      with open(filename) as f:
          for line in f:
              tmp = [x.strip() for x in line.split("=")]
              sym = tmp[0]
              val = [x.strip() for x in tmp[1].split()][0]
              setattr(self, sym, eval(val))

      """Keep a static copy of the parameters"""
      Params.instance = self

if __name__ == "__main__":
    P = Params("params.txt")
    for k in P.__dict__:
       print "{} : {}".format(k, P.__dict__[k])
