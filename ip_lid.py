from dolfin import *
from sys import argv


class NoSlipDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 1.0 - DOLFIN_EPS

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 1.0 - DOLFIN_EPS

def run_IP_lid( V_order="2",
                 h_num="8",
                 log_file="results/ip_lid.log",
                 ):
    V_element="CG",
    V_order = int(V_order);h_num = int(h_num);
    print V_order, h_num

    mesh = UnitSquare(h_num, h_num, "crossed");

    V = VectorFunctionSpace(mesh, "CG", V_order)

    noslip_domain = NoSlipDomain()
    noslip_val = Constant( (0.0, 0.0))
    top_domain = Top()
    top_val = Expression(("x[0]*(1.0 - x[0])", "0.0"))
    bc0 = DirichletBC(V, noslip_val, noslip_domain)
    bc1 = DirichletBC(V, top_val, top_domain)
    bc = [bc0, bc1]

    v = TestFunction(V)
    u = TrialFunction(V)
    f = Constant( (0.0, 0.0))

    w = Function(V)
    c = 1e3

    a = (inner(grad(v), grad(u)) - c*div(v)*div(u))*dx
    L = (inner(v, f) + inner(div(v), div(w)))*dx

    iters = 0
    max_iters = 100

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

    print "iters", iters, "U_m_u", U_m_u, "u_div", u_div
    dofs = len(U.vector())

    fp = open(log_file, 'a')
    fp.write("%(V_space)s, %(V_order)d,"\
             " %(h_num)d, %(dofs)d, %(time)e, %(iters)d, %(u_div)e\n" \
             % {"V_space":V_element,
                "V_order":V_order,
                "h_num":h_num,
                "dofs":dofs,
                "time":t.value(),
                "iters":iters,
                "u_div":u_div,
                })

def main():
    print argv
    run_IP_lid(*argv[1:])



if __name__ == "__main__":
    main()

