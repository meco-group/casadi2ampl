from casadi import *

def AMPLexport(nlp,data):

  fun = Function('nlp',nlp,["x","p"],["f","g"])
  fun = fun.expand()

  ins = fun.sx_in()
  fun = Function('f',ins,fun.call(ins),{"live_variables": False})

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
  
  p = DM(fun.sparsity_in(1))
  if "p" in data: p[:,:] = data["p"]

  main= ";\n".join(["var x%d := %.16f, >= %.16f, <= %.16f" % (i,float(x0[i]),float(lbx[i]),float(ubx[i])) for i in range(x0.shape[0])])+";\n"
  
  main+= ";\n".join(["param p%d := %.16f" % (i,float(p[i])) for i in range(p.shape[0])])+";\n"
  
  algorithm = "".join(str(fun).split("\n"))
  for i,l in enumerate(algorithm):
    if l.startswith("@"):
      break
  algorithm = "".join(algorithm[i:])
    
  algorithm = algorithm.replace("@","at")
  algorithm = re.sub("input\[0\]\[(\d+)\]",lambda m: "x%d" % (int(m.group(1))),algorithm)
  algorithm = re.sub("input\[1\]\[(\d+)\]",lambda m: "p%d" % (int(m.group(1))),algorithm)
  algorithm = re.sub("output\[0\]\[(\d+)\]",r"f",algorithm)
  algorithm = re.sub("output\[1\]\[(\d+)\]",r"g\1",algorithm)
  algorithm = re.sub(r"\bsq\((.*?)\)",r"(\1)^2",algorithm)
  
  is_var = dict()
  
  algorithms = []
  for a in algorithm.split(";")[:-1]:
    lhs, rhs = a.split("=")
    
    
    var = "x" in rhs

    if not var:
      for m in re.findall("(at\d+)",rhs):
        if is_var[m]:
          var = True
          break

    is_var[lhs.strip()] = var
    
    if var:
      algorithms.append("var " +a)
    else:
      algorithms.append("param " +a)
    
  main += ";\n".join(algorithms)+ ";\n"



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
  constr = ";\n".join(constr)+";\n"
      

  print """

reset;

{main}

minimize f1: f;

s.t.
{constr}

option solver knitro;

solve;

  """.format(main=main,constr=constr)
