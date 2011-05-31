from dolfin import *
from sys import argv

class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 1.0 - DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 1.0 - DOLFIN_EPS

class PinPoint(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS


def run_TH_lid( V_element="CG",
                 V_order="2",
                 P_element="CG",
                 P_order="1",
                 h_num="8",
                 log_file="results/mixed_lid.log",
                 stabilized=False
                 ):
    V_order = int(V_order); P_order = int(P_order); h_num = int(h_num);
    print V_order, P_order, h_num

    mesh = UnitSquare(h_num, h_num, "crossed");


    # Define function spaces
    if V_element.lower() == "mini":
        P = VectorFunctionSpace(mesh, "CG", V_order)
        B = VectorFunctionSpace(mesh, "Bubble", V_order + 2)
        V = P + B
    else:
        V = VectorFunctionSpace(mesh, V_element, V_order)
    Q = FunctionSpace(mesh, P_element, P_order)
    W = V * Q
    VP10 = VectorFunctionSpace(mesh, "CG", 10)
    P10 = FunctionSpace(mesh, "CG", 10)

    # No-slip boundary condition for velocity
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
    if stabilized:
        h = CellSize(mesh)
        beta = 0.2
        delta = beta*h*h
        a = (inner(grad(v), grad(u)) - div(v)*p + q*div(u) +\
             delta*inner(grad(q), grad(p)))*dx
        L = inner(v + delta*grad(q), f)*dx
    else:
        a = (inner(grad(v), grad(u)) - div(v)*p + q*div(u))*dx
        L = inner(v, f)*dx

    # Compute solution
    t = Timer("Solve timing");
    t.start();
    pde = VariationalProblem(a, L, bc)
    U = pde.solve()
    t.stop();
    dofs = len(U.vector())
    # Split the mixed solution using deepcopy
    # (needed for further computation on coefficient vector)
    (u, p) = U.split(True)

    Mdiv = div(u)*div(u)*dx
    v_div = assemble(Mdiv, mesh=mesh)
    print "Error of div: %.15g" % v_div

    fp = open(log_file, 'a')
    fp.write("%(V_space)s, %(V_order)d, %(P_space)s,"\
             " %(P_order)d, %(h_num)d, %(dofs)d, %(time)e, %(v_div)e\n" \
             % {"V_space":V_element,
                "V_order":V_order,
                "P_space":P_element,
                "P_order":P_order,
                "h_num":h_num,
                "dofs":dofs,
                "time":t.value(),
                "v_div":v_div,
                })
#     ufile_pvd = File("results/velocity_lid.pvd")
#     ufile_pvd << u
#     pfile_pvd = File("results/pressure_lid.pvd")
#     pfile_pvd << p

def main():
    print argv
    run_TH_lid(*argv[1:])


if __name__ == "__main__":
    main()

