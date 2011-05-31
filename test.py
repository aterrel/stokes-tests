from dolfin import *

class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

mesh = UnitSquare(8, 8, "crossed");

V = VectorFunctionSpace(mesh, "CG", 2)

noslip_domain = NoSlipDomain()
noslip = Expression(("x[0] * x[0] * x[1]", "-x[0] * x[1] * x[1]"), V=V)
bc = DirichletBC(V, noslip, noslip_domain)
    
v = TestFunction(V)
u = TrialFunction(V)
f = Expression(("-2*x[1] + 1.0", "2*x[0] + 1.0"), V=V)

w = Function(V)
c = 100

a = (inner(grad(v), grad(u)) - c*div(v)*div(u))*dx
L = (inner(v, f) + inner(div(v), div(w)))*dx

iters = 0
max_iters = 5

problem = VariationalProblem(a, L, bc)
U = problem.solve()
M = div(U)*dx
div_u = sqrt(assemble(M, mesh=mesh))

while iters < max_iters and div_u > DOLFIN_SQRT_EPS:
    w.vector().axpy(c, U.vector())
    U = problem.solve()
    div_u = sqrt(assemble(M, mesh=mesh))
    iters += 1
    print "div_u", div_u

print "iters", iters
print "div_u", sqrt(assemble(M, mesh=mesh))

