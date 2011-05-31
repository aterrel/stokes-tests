from dolfin import *

class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 1.0 - DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 1.0 - DOLFIN_EPS

mesh = UnitSquare(32, 32, "crossed");

# Define function spaces
V = VectorFunctionSpace(mesh, "Crouzeix-Raviart", 1)
Q = FunctionSpace(mesh, "DG", 0)
W = V + Q

noslip_domain = NoSlipDomain()
noslip = Constant(mesh, (0,0))
bc0 = DirichletBC(V, noslip, noslip_domain)
top = Top()
top_fun = Expression(("x[0]*(1.0 - x[0])","0.0"), V=V)
bc1 = DirichletBC(V, top_fun, top)
bc = [bc0, bc1]
    
(v, q) = TestFunctions(W)
du = TrialFunction(W)
u_m = Function(W)
u, p = split(u_m)
w = Function(V)
f = Constant(mesh, (0.0, 0.0))

L = (inner(grad(v), grad(u))  + inner(v, dot(w, grad(u))) - div(v)*p + q*div(u))*dx + inner(v, f)*dx

a = derivative(L, u_m, du)
    
# Compute solution
t = Timer("Solve timing");
t.start();
pde = VariationalProblem(a, L, bc, nonlinear=True)
U = pde.solve()
t.stop();
dofs = len(U.vector())
# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u, p) = U.split(True)

print "Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2")
print "Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2")
Mdiv = div(u)*div(u)*dx
v_div = assemble(Mdiv, mesh=mesh)
print "Error of div: %.15g" % v_div
