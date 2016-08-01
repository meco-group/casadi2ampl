from casadi import *
from AMPLexport import AMPLexport

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

AMPLexport(nlp,data)
