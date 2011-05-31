from dolfin import *
from sys import argv

problems = [
    {'f' : ("-2*x[1] + 1.0", "2*x[0] + 1.0"),
     'u' : ("x[0] * x[0] * x[1]", "-x[0] * x[1] * x[1]"),
     'p' : "x[0] + x[1] - 1.0"
     },
    {'f' : ("6*pi*pi*2*sin(2*pi*x[0])*cos(2*pi*x[1])",
            "-10*pi*pi*cos(2*pi*x[0])*sin(2*pi*x[1])"),
     'u' : ("sin(2*pi*x[0])*cos(2*pi*x[1])",
            "-cos(2*pi*x[0])*sin(2*pi*x[1])"),
     'p' : "pi*cos(2*pi*x[0])*cos(2*pi*x[1])"
     },
    {'f' : ("28*pi*pi*sin(4*pi*x[0])*cos(4*pi*x[1])",
            "-36*pi*pi*cos(4*pi*x[0])*sin(4*pi*x[1])"),
     'u' : ("sin(4*pi*x[0])*cos(4*pi*x[1])",
            "-cos(4*pi*x[0])*sin(4*pi*x[1])"),
     'p' : "pi*cos(4*pi*x[0])*cos(4*pi*x[1])"
     }
    ]

# Subdomain for boundary
class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class PinPoint(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS


def run_TH_test( V_element="CG",
                 V_order="2",
                 P_element="CG",
                 P_order="1",
                 h_num="8",
                 problem="0",
                 log_file="results/mixed.log",
                 stabilized=False
                 ):
    V_order = int(V_order); P_order = int(P_order); h_num = int(h_num);
    prob = problems[int(problem)]
    print V_order, P_order, h_num, prob

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
    noslip = Expression(prob['u'])
    pinpoint = PinPoint()
    pin_val = Expression(prob['p'])

    bc0 = DirichletBC(W.sub(0), noslip, noslip_domain)
    bc1 = DirichletBC(W.sub(1), pin_val, pinpoint, "pointwise")
    bc = [bc0, bc1]

    # Define variational problem
    (v, q) = TestFunctions(W)
    (u, p) = TrialFunctions(W)
    f = Expression(prob['f'])
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

    print "Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2")
    print "Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2")

    u_ex = Expression(prob['u'], element=VP10.ufl_element())
    M = inner((u_ex - u),(u_ex - u))*dx
    v_err = assemble(M, mesh=mesh)
    print "Error of velocity coefficient vector: %.15g" % v_err

    p_ex = Expression(prob['p'], element=P10.ufl_element())
    Mp = (p_ex - p)*(p_ex - p)*dx
    p_err = assemble(Mp, mesh=mesh)
    print "Error of pressure coefficient vector: %.15g" % p_err

    Mdiv = div(u)*div(u)*dx
    v_div = assemble(Mdiv, mesh=mesh)
    print "Error of div: %.15g" % v_div

    fp = open(log_file, 'a')
    fp.write("%(problem)d, %(V_space)s, %(V_order)d, %(P_space)s,"\
             " %(P_order)d, %(h_num)d, %(dofs)d, %(time)e, %(v_err)e, %(p_err)e, %(v_div)e\n" \
             % {"problem":int(problem),
                "V_space":V_element,
                "V_order":V_order,
                "P_space":P_element,
                "P_order":P_order,
                "h_num":h_num,
                "dofs":dofs,
                "time":t.value(),
                "v_err":v_err,
                "p_err":p_err,
                "v_div":v_div,
                })


def main():
    print argv
    run_TH_test(*argv[1:])


if __name__ == "__main__":
    main()

