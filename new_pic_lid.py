from dolfin import *

class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 1.0 - DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 1.0 - DOLFIN_EPS

class PinPoint(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS

V_order=2
P_order=1
h_num=16
mesh = UnitSquare(h_num, h_num, "crossed");
V_element="CG"
V = VectorFunctionSpace(mesh, V_element, V_order)
P_element="CG"
Q = FunctionSpace(mesh, P_element, P_order)
W = V*Q
noslip_domain = NoSlipDomain()
noslip_val = Constant((0.0, 0.0))
top_domain = Top()
top_val = Expression(("x[0]*(1.0 - x[0])", "0.0"))
pinpoint = PinPoint()
pin_val = Constant(0.0)
bc0 = DirichletBC(W.sub(0), noslip_val, noslip_domain)
bc1 = DirichletBC(W.sub(0), top_val, top_domain)
bc2 = DirichletBC(W.sub(1), pin_val, pinpoint, "pointwise")
bc = [bc0, bc1, bc2]

# Define variational problem
(v, q) = TestFunctions(W)
(u, p) = TrialFunctions(W)
f = Constant((0.0, 0.0))
a = (inner(grad(v), grad(u)) - div(v)*p + q*div(u))*dx
L = inner(v, f)*dx

# Compute solution
pde = VariationalProblem(a, L, bc)
U = pde.solve()

# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u0, p0) = U.split(True)

a_div = inner(p, q)*dx
L_div = inner(q, div(u0))
div_u = VariationalProblem(a_div, L_div).solve()

file = File("TH_div.pvd")
file << div_u


