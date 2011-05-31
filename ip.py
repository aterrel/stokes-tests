from dolfin import *
from sys import argv
from mixed import problems, NoSlipDomain


def run_IP_test( V_order="2",
                 h_num="8",
                 problem="0",
                 log_file="results/ip_test.log",
                 ):
    V_element="CG",
    V_order = int(V_order);h_num = int(h_num);
    prob = problems[int(problem)]
    print V_order, h_num, prob

    mesh = UnitSquare(h_num, h_num, "crossed");

    V = VectorFunctionSpace(mesh, "CG", V_order)
    VP10 = VectorFunctionSpace(mesh, "CG", 10)
    P10 = FunctionSpace(mesh, "CG", 10)

    noslip_domain = NoSlipDomain()
    noslip = Expression(prob['u'])
    bc = DirichletBC(V, noslip, noslip_domain)

    v = TestFunction(V)
    u = TrialFunction(V)
    f = Expression(prob['f'])

    w = Function(V)
    c = 100

    a = (inner(grad(v), grad(u)) - c*div(v)*div(u))*dx
    L = (inner(v, f) + inner(div(v), div(w)))*dx

    iters = 0
    max_iters = 15

    pde = VariationalProblem(a, L, bc)
    t = Timer("timer")
    t.start()
    U_m_u = 1

    while iters < max_iters and U_m_u > 1e-8: # and u_div > DOLFIN_EPS:
        U = pde.solve()
        w.vector().axpy(c, U.vector())
        if iters == 0:
            U_m_u = 1
        else:
            U_m_u = (U.vector() - u_old_vec).norm('l2')
        u_old_vec = U.vector().copy()
        iters += 1


    t.stop()
    M = div(U)*div(U)*dx
    u_div = assemble(M, mesh=mesh)

    print "iters", iters
    print "U_m_u", U_m_u
    print "u_div", u_div
    dofs = len(U.vector())

    u_ex = Expression(prob['u'], element=VP10.ufl_element())
    M = inner((u_ex - U),(u_ex - U))*dx
    v_err = assemble(M, mesh=mesh)
    print "Error of velocity coefficient vector: %.15g" % v_err

    W = FunctionSpace(mesh, "CG", V_order-1)
    q = TestFunction(W)
    p = TrialFunction(W)
    pde_1 = VariationalProblem(inner(q,p)*dx, inner(q,div(w))*dx)
    P = pde_1.solve()
    ave_p = P.vector().sum()/ len(P.vector())
    for i, v in enumerate(P.vector()):
        P.vector()[i] = v - ave_p

    p_ex = Expression(prob['p'], element=P10.ufl_element())
    Mp = (p_ex - P)*(p_ex - P)*dx
    p_err = assemble(Mp, mesh=mesh)
    print "Error of pressure coefficient vector: %.15g" % p_err

    fp = open(log_file, 'a')
    fp.write("%(problem)d, %(V_space)s, %(V_order)d,"\
             " %(h_num)d, %(dofs)d, %(time)e, %(iters)d, %(v_err)e, %(p_err)e, %(u_div)e\n" \
             % {"problem":int(problem),
                "V_space":V_element,
                "V_order":V_order,
                "h_num":h_num,
                "dofs":dofs,
                "time":t.value(),
                "iters":iters,
                "v_err":v_err,
                "p_err":p_err,
                "u_div":u_div,
                })

def main():
    print argv
    run_IP_test(*argv[1:])



if __name__ == "__main__":
    main()

