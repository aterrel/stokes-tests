from dolfin import *
from sys import argv

prob = {'f' : ("28*pi*pi*sin(4*pi*x[0])*cos(4*pi*x[1])",
               "-36*pi*pi*cos(4*pi*x[0])*sin(4*pi*x[1])"),
        'u' : ("sin(4*pi*x[0])*cos(4*pi*x[1])",
               "-cos(4*pi*x[0])*sin(4*pi*x[1])"),
        'p' : "pi*cos(4*pi*x[0])*cos(4*pi*x[1])"
        }


# Subdomain for boundary
class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class PinPoint(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS


def run_mix_test(V_order=2, h_num=8):
    mesh = UnitSquare(h_num, h_num, "crossed");

    # Define function spaces
    P = VectorFunctionSpace(mesh, "CG", V_order)
    B = VectorFunctionSpace(mesh, "Bubble", V_order + 2)
    V = P + B

    Q = FunctionSpace(mesh, "CG", V_order)
    R = FunctionSpace(mesh, "R", 0)

    W = MixedFunctionSpace([V, Q, R])

    VP10 = VectorFunctionSpace(mesh, "CG", 10)
    P10 = FunctionSpace(mesh, "CG", 10)

    # No-slip boundary condition for velocity
    noslip_domain = NoSlipDomain()
    noslip = Expression(prob['u'])
    pinpoint = PinPoint()
    pin_val = Expression(prob['p'])

    bc0 = DirichletBC(W.sub(0), noslip, noslip_domain)
    bc1 = DirichletBC(W.sub(1), pin_val, pinpoint, "pointwise")
    bc = [bc0, bc1]

    # Define variational problem
    (v, q, d) = TestFunctions(W)
    (u, p, c) = TrialFunctions(W)
    f = Expression(prob['f'])
    a = (inner(grad(v), grad(u)) - div(v)*p + q*div(u) + q*c + d*p)*dx
    L = inner(v, f)*dx

    # Compute solution
    pde = VariationalProblem(a, L, bc)
    U = pde.solve()
    # Split the mixed solution using deepcopy
    # (needed for further computation on coefficient vector)
    (u, p, _) = U.split(True)

    u_ex = Expression(prob['u'], element=VP10.ufl_element())
    M = inner((u_ex - u),(u_ex - u))*dx
    v_err = assemble(M, mesh=mesh)

    p_ex = Expression(prob['p'], element=P10.ufl_element())
    Mp = (p_ex - p)*(p_ex - p)*dx
    p_err = assemble(Mp, mesh=mesh)

    Mdiv = div(u)*div(u)*dx
    div_err = assemble(Mdiv, mesh=mesh)

    return v_err, p_err, div_err

def main():
    errs = {}
    for order in xrange(1, 4):
        N = 4
        errs[order] = {}
        while N < 18:
            errs[order][N] = run_mix_test(V_order=order, h_num=N)
            N *= 2

    print "Velocity order, N, v_err, p_err, div_err"
    for order in xrange(1, 4):
        N = 4
        while N < 18:
            print "%d, %d, %.2e, %.2e, %.2e" % \
                  (order, N, errs[order][N][0], errs[order][N][1], errs[order][N][2])
            N *= 2


if __name__ == "__main__":
    main()

