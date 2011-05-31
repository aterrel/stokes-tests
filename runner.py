from sys import argv
from mixed import run_TH_test
from ip import run_IP_test
from mixed_lid import run_TH_lid
from ip_lid import run_IP_lid


HMIN = 4
HMAX = 64
HCHEAPMAX = 128

print len(argv), "Args:", argv

if len(argv) == 1:
    argv.append("ALL")

if argv[1] in ["ALL", "TH"]:
    print "Running TH"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(2,4):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_order = i,
                        P_order = i-1,
                        log_file = "results/TH.log")
        h_mesh *= 2

if argv[1] in ["ALL", "TH_STAB"]:
    print "Running TH_STAB"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(2,6):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_order = i,
                        P_order = i-1,
                        log_file = "results/TH_STAB.log",
                        stabilized=True)
        h_mesh *= 2


if argv[1] in ["ALL", "CR"]:
    print "Running CR"
    h_mesh = HMIN
    while h_mesh <= HCHEAPMAX:
        run_TH_test(problem = 2,
                    h_num = h_mesh,
                    V_element="Crouzeix-Raviart",
                    V_order = 1,
                    P_element = "DG",
                    P_order = 0,
                    log_file = "results/CR.log")
        h_mesh *= 2

# if argv[1] in ["ALL", "CPCP"]:
#     print "Running CPCP"
#     h_mesh = HMIN
#     while h_mesh <= HMAX:
#         for i in xrange(2,6):
#             run_TH_test(problem = 2,
#                         h_num = h_mesh,
#                         V_order = i,
#                         P_element = "DG",
#                         P_order = i-1,
#                         log_file = "results/CPCP.log")
#         h_mesh *= 2

if argv[1] in ["ALL", "BF"]:
    print "Running BF"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(2,6):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_order = i,
                        P_element = "DG",
                        P_order = i-2,
                        log_file = "results/BF.log")
        h_mesh *= 2

if argv[1] in ["ALL", "MINI"]:
    print "Running MINI"
    h_mesh = HMIN
    while h_mesh <= HCHEAPMAX:
        for i in xrange(1,4):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_element = "MINI",
                        V_order = i,
                        P_order = i,
                        log_file = "results/MINI.log")
        h_mesh *= 2

if argv[1] in ["ALL", "MINISTAB"]:
    print "Running MINISTAB"
    h_mesh = HMIN
    while h_mesh <= HCHEAPMAX:
        for i in xrange(1,4):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_element = "MINI",
                        V_order = i,
                        P_order = i,
                        log_file = "results/MINIStab.log",
                        stabilized = True)
        h_mesh *= 2


if argv[1] in ["ALL", "STAB"]:
    print "Running STAB"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(1,6):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_order = i,
                        P_element = "CG",
                        P_order = i,
                        log_file = "results/STAB.log",
                        stabilized = True)
        h_mesh *= 2

if argv[1] in ["ALL", "IP"]:
    print "Running IP"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(4,6):
            run_IP_test(problem = 2,
                        h_num = h_mesh,
                        V_order = i,
                        log_file = "results/IP.log")
        h_mesh *= 2

if argv[1] in ["ALL", "TH_lid"]:
    print "Running TH_lid"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(2,6):
            run_TH_lid(h_num = h_mesh,
                       V_order = i,
                       P_order = i-1,
                       log_file = "results/TH_lid.log")
        h_mesh *= 2

if argv[1] in ["ALL", "TH_STAB_lid"]:
    print "Running TH_STAB_lid"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(2,6):
            run_TH_test(problem = 2,
                        h_num = h_mesh,
                        V_order = i,
                        P_order = i-1,
                        log_file = "results/TH_STAB_lid.log",
                        stabilized=True)
        h_mesh *= 2


if argv[1] in ["ALL", "CR_lid"]:
    print "Running CR_lid"
    h_mesh = HMIN
    while h_mesh <= HCHEAPMAX:
        run_TH_lid(h_num = h_mesh,
                   V_element="Crouzeix-Raviart",
                   V_order = 1,
                   P_element = "DG",
                   P_order = 0,
                   log_file = "results/CR_lid.log")
        h_mesh *= 2

# if argv[1] in ["ALL", "CPCP_lid"]:
#     print "Running CPCP_lid"
#     h_mesh = HMIN
#     while h_mesh <= HMAX:
#         for i in xrange(2,6):
#             run_TH_lid(h_num = h_mesh,
#                        V_order = i,
#                        P_element = "DG",
#                        P_order = i-1,
#                        log_file = "results/CPCP_lid.log")
#         h_mesh *= 2

if argv[1] in ["ALL", "BF_lid"]:
    print "Running BF_lid"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(2,6):
            run_TH_lid(h_num = h_mesh,
                       V_order = i,
                       P_element = "DG",
                       P_order = i-2,
                       log_file = "results/BF_lid.log")
        h_mesh *= 2


if argv[1] in ["ALL", "MINISTAB_lid"]:
    print "Running MINISTAB_lid"
    h_mesh = HMIN
    while h_mesh <= HCHEAPMAX:
        for i in xrange(1,4):
            run_TH_lid(h_num = h_mesh,
                       V_element = "MINI",
                       V_order = i,
                       P_order = i,
                       log_file = "results/MINIStab_lid.log",
                       stabilized = True)
        h_mesh *= 2

if argv[1] in ["ALL", "MINI_lid"]:
    print "Running MINI_lid"
    h_mesh = HMIN
    while h_mesh <= HCHEAPMAX:
        for i in xrange(1,4):
            run_TH_lid(h_num = h_mesh,
                       V_element = "MINI",
                       V_order = i,
                       P_order = i,
                       log_file = "results/MINI_lid.log")
        h_mesh *= 2


if argv[1] in ["ALL", "STAB_lid"]:
    print "Running STAB_lid"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(1,6):
            run_TH_lid(h_num = h_mesh,
                       V_order = i,
                       P_element = "CG",
                       P_order = i,
                       log_file = "results/STAB_lid.log",
                       stabilized = True)
        h_mesh *= 2

if argv[1] in ["ALL", "IP_lid"]:
    print "Running IP_lid"
    h_mesh = HMIN
    while h_mesh <= HMAX:
        for i in xrange(4,6):
            run_IP_lid(h_num = h_mesh,
                       V_order = i,
                       log_file = "results/IP_lid.log")
        h_mesh *= 2

