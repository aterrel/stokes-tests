from dolfin import *

class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 1.0 - DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 1.0 - DOLFIN_EPS

mesh = UnitSquare(32, 32, "crossed");

V = VectorFunctionSpace(mesh, "CG", 2)

noslip_domain = NoSlipDomain()
noslip = Expression(("x[0] * x[0] * x[1]", "-x[0] * x[1] * x[1]"), V=V)
bc0 = DirichletBC(V, noslip, noslip_domain)
top = Top()
top_fun = Expression(("x[0]*(1.0 - x[0])","0.0"), V=V)
bc1 = DirichletBC(V, top_fun, top)
bc = [bc0, bc1]
    
v = TestFunction(V)
u = TrialFunction(V)
f = Constant(mesh, (0.0, 0.0))

w = Function(V)
c = 1e8

a = (inner(grad(v), grad(u)) - c*div(v)*div(u))*dx
L = (inner(v, f) + inner(div(v), div(w)))*dx

iters = 0
max_iters = 5

problem = VariationalProblem(a, L, bc)
U = problem.solve()
w.vector().axpy(c, U.vector())
M = div(U)*div(U)*dx
div_u = assemble(M, mesh=mesh)

while iters < max_iters and div_u > DOLFIN_EPS:
    print w.vector().norm('l2'), U.vector().norm('l2')
    U = problem.solve()
    w.vector().axpy(c, U.vector())
    div_u = assemble(M, mesh=mesh)
    iters += 1
    print "div_u", div_u

print "iters", iters
print "div_u", assemble(M, mesh=mesh)

