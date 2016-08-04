from casadi import *
import numpy as np

def AMPLexport(nlp,data):
  def isinf(a):
    return bool(a>=1e20) or bool(a<=-1e20)

  fun = Function('nlp',nlp,["x","p"],["f","g"])
  fun = fun.expand()

  ins = fun.sx_in()
  fun = Function('f',ins,fun.call(ins),{"live_variables": False})

  lbg = -np.inf*DM(fun.sparsity_out(1))
  ubg = np.inf*DM(fun.sparsity_out(1))
  if "lbg" in data: lbg[:,:] = data["lbg"]
  if "ubg" in data: ubg[:,:] = data["ubg"]


  lbx = -np.inf*DM(fun.sparsity_in(0))
  x0 = DM(fun.sparsity_in(0))
  ubx = np.inf*DM(fun.sparsity_in(0))
  if "lbx" in data: lbx[:,:] = data["lbx"]
  if "ubx" in data: ubx[:,:] = data["ubx"]
  if "x0" in data: x0[:,:] = data["x0"]
  
  p = DM(fun.sparsity_in(1))
  if "p" in data: p[:,:] = data["p"]
  main = ""
  for i in range(p.shape[0]):
     main += "param p%d = %.16f;\n" % (i,float(p[i]))
  for i in range(x0.shape[0]):
    if isinf(lbx[i]) and isinf(ubx[i]):
      main+="var x%d := %.16f;\n" % (i,float(x0[i]))
    elif isinf(lbx[i]):
      main+="var x%d := %.16f, <= %.16f;\n" % (i,float(x0[i]),float(ubx[i]))
    elif isinf(ubx[i]):
      main+="var x%d := %.16f, >= %.16f;\n" % (i,float(x0[i]),float(lbx[i]))
    else:
      main+="var x%d := %.16f, >= %.16f, <= %.16f;\n" % (i,float(x0[i]),float(lbx[i]),float(ubx[i]))
  
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
    try:
      lhs, rhs = a.split("=",1)
    except:
      print 'There was an error here'
    
    print rhs
    if ">" in rhs or "<" in rhs or "=" in rhs:
      rhs = "if %s then 1 else 0" % rhs

    var = "x" in rhs

    if not var:
      for m in re.findall("(at\d+)",rhs):
        if is_var[m]:
          var = True
          break

    is_var[lhs.strip()] = var
    
    a = "%s = %s" % (lhs, rhs)
    if var:
      algorithms.append("var " +a)
    else:
      algorithms.append("param " +a)
    
  main += ";\n".join(algorithms)+ ";\n"



  constr = []


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
      constr.append("  con%d_upper: g%i >= %.16f" % (i,i,float(ubg[i])))
  constr = ";\n".join(constr)+";\n"
  
  displayvar = ",".join("x%d" % (i) for i in range(x0.shape[0]))

  print """

reset;

{main}

minimize f1: f;

s.t.
{constr}

option solver knitro;

solve;

display {displayvar};

  """.format(main=main,constr=constr,displayvar=displayvar)
