from dolfin import *
from sys import argv

problems = [
    {'f' : ("-2*x[1] + 1.0", "2*x[0] + 1.0"),
     'u' : ("x[0] * x[0] * x[1]", "-x[0] * x[1] * x[1]"),
     'p' : "x[0] + x[1] - 1.0"
     },
    {'f' : ("8*pi*pi*2*sin(2*pi*x[0])*cos(2*pi*x[1])",
            "-8*pi*pi*cos(2*pi*x[0])*sin(2*pi*x[1])"),
     'u' : ("sin(2*pi*x[0])*cos(2*pi*x[1])",
            "-cos(2*pi*x[0])*sin(2*pi*x[1])"),
     'p' : "sin(2*pi*x[0])*sin(2*pi*x[0])"
     },
    {'f' : ("18*pi*pi*sin(3*pi*x[0])*cos(3*pi*x[1])",
            "-18*pi*pi*cos(3*pi*x[0])*sin(3*pi*x[1])"),
     'u' : ("sin(3*pi*x[0])*cos(3*pi*x[1])",
            "-cos(3*pi*x[0])*sin(3*pi*x[1])"),
     'p' : "sin(3*pi*x[0])*sin(3*pi*x[1])"
     }
    ]

# Subdomain for boundary
class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

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
    class PinPoint(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0] < 1./h_num + DOLFIN_EPS and \
                   x[1] < 1./h_num + DOLFIN_EPS

    mesh = UnitSquare(h_num, h_num, "crossed");

    sub_domains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)
    sub_domains.set_all(3)
    noslip_domain = NoSlipDomain()
    noslip_domain.mark(sub_domains, 0)
    pinpoint = PinPoint()
    pinpoint.mark(sub_domains, 1)

    # Define function spaces
    V = FunctionSpace(mesh, V_element, V_order)
    Q = FunctionSpace(mesh, P_element, P_order)
    W = V + Q
    VP10 = VectorFunctionSpace(mesh, "CG", 10)
    P10 = FunctionSpace(mesh, "CG", 10)

    # No-slip boundary condition for velocity
    noslip = Expression(prob['u'], V=V)
    bc0 = DirichletBC(W.sub(0), noslip, sub_domains, 0)
    
    # Boundary condition for pressure at pinpoint
    #zero = Constant(mesh, 0.0)
    #bc1 = DirichletBC(W.sub(1), zero, sub_domains, 1), "pointwise")
    pin_val = Expression(prob['p'], V=Q)
    bc1 = DirichletBC(W.sub(1), pin_val, sub_domains, 1)

    # Collect boundary conditions
    bc = [bc0, bc1]

    # Define variational problem
    (v, q) = TestFunctions(W)
    (u, p) = TrialFunctions(W)
    f = Expression(prob['f'], V=V)
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
    pde = VariationalProblem(a, L, bc)
    U = pde.solve()
    
    # Split the mixed solution using deepcopy
    # (needed for further computation on coefficient vector)
    (u, p) = U.split(True)
    
    print "Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2")
    print "Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2")
    

    u_ex = Expression(prob['u'], V=VP10)
    M = inner((u_ex - u),(u_ex - u))*dx
    v_err = assemble(M, mesh=mesh)
    #VPu = project(u_ex - u, VP10)
    #print "Error of velocity coefficient vector: %.15g" % VPu.vector().norm("l2")
    print "Error of velocity coefficient vector: %.15g" % v_err


    p_ex = Expression(prob['p'], V=P10)
    Mp = (p_ex - p)*(p_ex - p)*dx
    p_err = assemble(Mp, mesh=mesh)
#    Pp = project(p_ex - p, P10)
#    print "Error of pressure coefficient vector: %.15g" % Pp.vector().norm("l2")
    print "Error of pressure coefficient vector: %.15g" % p_err

    Mdiv = div(u)*div(u)*dx
    v_div = assemble(Mdiv, mesh=mesh)

    fp = open(log_file, 'a')
    fp.write("%(problem)d, %(V_space)s, %(V_order)d, %(P_space)s,"\
             " %(P_order)d, %(h_num)d, %(v_err)e, %(p_err)e, %(v_div)e\n" \
             % {"problem":int(problem),
                "V_space":V_element,
                "V_order":V_order,
                "P_space":P_element,
                "P_order":P_order,
                "h_num":h_num,
                "v_err":v_err,
                "p_err":p_err,
                "v_div":v_div,
                })
    
def main():
    print argv        
    run_TH_test(*argv[1:])



if __name__ == "__main__":
    main()

