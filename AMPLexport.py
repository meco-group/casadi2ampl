from casadi import *


# Declare variables
x = SX.sym("x")
y = SX.sym("y")
z = SX.sym("z")

# Formulate the NLP
f = x**2 + 100*z**2
g = z + (1-x)**2 - y
nlp = {'x':vertcat(x,y,z), 'f':f, 'g':g}

# Create an NLP solver
solver = nlpsol("solver", "ipopt", nlp)

data = {"x0":[2.5,3.0,0.75],"lbg":0,"ubg":0}

def AMPLexport(solver,data):

  fun = solver.oracle()
  fun = fun.expand()

  ins = fun.sx_in()
  fun = Function('f',ins,fun.call(ins),{"live_variables": False})

  main = """
var x{1..%d} := 0;
var p{1..%d} := 0;
  """ % (fun.nnz_in(0),fun.nnz_in(1))


  algorithm = "".join(str(fun).split("\n"))
  for i,l in enumerate(algorithm):
    if l.startswith("@"):
      break
  algorithm = "".join(algorithm[i:])
    
  algorithm = algorithm.replace("@","at")
  algorithm = re.sub("input\[0\]\[(\d+)\]",lambda m: "x[%d]" % (int(m.group(1))+1),algorithm)
  algorithm = re.sub("output\[0\]\[(\d+)\]",r"f",algorithm)
  algorithm = re.sub("output\[1\]\[(\d+)\]",r"g\1",algorithm)
  main += ";\n".join([ "var "+s for s in  algorithm.split(";")[:-1]])

  lbg = DM(fun.sparsity_out(1))
  ubg = DM(fun.sparsity_out(1))
  if "lbg" in data: lbg[:,:] = data["lbg"]
  if "ubg" in data: ubg[:,:] = data["ubg"]

  lbx = DM(fun.sparsity_in(0))
  x0 = DM(fun.sparsity_in(0))
  ubx = DM(fun.sparsity_in(0))
  if "lbx" in data: lbx[:,:] = data["lbx"]
  if "ubx" in data: ubx[:,:] = data["ubx"]
  if "x0" in data: x0[:,:] = data["x0"]


  bounds = ";\n".join(["let x[%d] := %.16f, >= %.16f, <= %.16f" % (i,float(x0[i]),float(lbx[i]),float(ubx[i])) for i in range(x0.shape[0])])

  constr = []

  def isinf(a):
    return bool(a>=1e-20) or bool(a<=-1e-20)

  for i in range(lbg.shape[0]):
    if isinf(lbg[i]) and isinf(ubg[i]):
      continue
    elif lbg[i]==ubg[i]:
      constr.append("  con%d: g%i = %.16f" % (i,i,float(lbg[i])))
    elif isinf(lbg[i]):
      constr.append("  con%d: g%i <= %.16f" % (i,i,float(ubg[i])))
    elif isinf(ubg[i]):
      constr.append("  con%d: g%i >= %.16f" % (i,i,float(lbg[i])))
    else:
      constr.append("  con%d_lower: g%i >= %.16f" % (i,i,float(lbg[i])))
      constr.append("  con%d_upper: g%i >= %.16f" % (i,i,float(lbg[i])))
  constr = ";\n".join(constr)
      

  print """

reset;

{main}
{bounds}

minimize f: f;

s.t.
{constr}

option solver knitro;

solve;

display x;

  """.format(main=main,constr=constr,bounds=bounds)
