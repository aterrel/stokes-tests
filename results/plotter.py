import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

#mixed_files = ["TH", "TH_STAB", "STAB", "CD", "MINI", "MINI_STAB", "CR"]
mixed_files = ["TH", "STAB", "CD", "MINI", "CR"]
ip_files = ["SV"]

def read_mixed_log(filename, test="test"):
#    print "Reading", filename
    fp = open(filename, 'r')
    vals = []
    for line in fp.readlines():
        line = line.strip()
        if not line:
            continue
        if test=="ldc":
            _, order, _, _, n, _, rt, div = line.split(",")
            vals.append({"order":int(order), "n":int(n), "rt":float(rt),
                         "div":float(div)})
        elif test=="test":
            _, _, order, _, _, n, _, rt, v, p, div = line.split(",")
            vals.append({"order":int(order), "n":int(n), "rt":float(rt),
                         "div":float(div), "v":float(v), "p":float(p)})
    return vals

def read_ip_log(filename, test="test"):
#    print "Reading", filename
    fp = open(filename, 'r')
    vals = []
    for line in fp.readlines():
        line = line.strip()
        if not line:
            continue
        if test=="ldc":
            _, _, order, n, _, rt, _, div = line.split(",")
            vals.append({"order":int(order), "n":int(n), "rt":float(rt),
                         "div":float(div)})
        elif test=="test":
            _, _, _, order, n, _, rt, _, v, p, div = line.split(",")
            vals.append({"order":int(order), "n":int(n), "rt":float(rt),
                         "div":float(div), "v":float(v), "p":float(p)})
    return vals

def create_4th_order_plots():
    vel_dict = {}
    pressure_dict = {}
    div_dict = {}
    runtime_dict = {}
    meshes = [4, 8, 16, 32, 64]
    for f in mixed_files + ip_files:
        print "Processing", f
        if f in ip_files:
            vals = read_ip_log(f+".log")
        else:
            vals = read_mixed_log(f+".log")
#        print "len(vals)", len(vals)
#        print vals
        if f in ["MINI", "CR", "MINI_STAB"]:
            scaled_meshes = [2*m for m in meshes]
            vals = filter(lambda v: v["order"] == 1 and v["n"] in scaled_meshes, vals)
            uniq_vals = []
            for n in scaled_meshes:
                for v in vals:
                    if v["n"] == n:
                        uniq_vals.append(v)
                        break
        else:
            vals = filter(lambda v: v["order"] == 4 and v["n"] in meshes, vals)
            uniq_vals = []
            for n in meshes:
                for v in vals:
                    if v["n"] == n:
                        uniq_vals.append(v)
                        break
        print "len(uniq_vals)", len(uniq_vals)
        vel_dict[f] = [v["v"] for v in uniq_vals]
        pressure_dict[f] = [v["p"] for v in uniq_vals]
        div_dict[f] = [v["div"] for v in uniq_vals]
        runtime_dict[f] = [v["rt"] for v in uniq_vals]
    bargraph_from_dict(vel_dict, meshes, log=True, save="vel_4.pdf",
                       title="4th Order Velocity Error", xlabel="Mesh size(n)",
                       ylabel="Velocity Error ($L^2$)", legend=True, loc='lower left')
    bargraph_from_dict(pressure_dict, meshes, log=True, save="press_4.pdf",
                       title="4th Order Pressure Error", xlabel="Mesh size(n)",
                       ylabel="Pressure Error ($L^2$)", legend=True, loc='lower left')
    bargraph_from_dict(div_dict, meshes, log=True, save="div_4.pdf",
                       title="4th Order Divergence Error", xlabel="Mesh size(n)",
                       ylabel="Divergence Error ($L^2$)", legend=True, loc='center left')
    bargraph_from_dict(runtime_dict, meshes, log=True, save="run_4.pdf",
                       title="4th Order Runtime", xlabel="Mesh size(n)",
                       ylabel="Runtime (s)", legend=True)


    for f in mixed_files+ip_files:
        print "Processing", f+"_lid"
        if f in ip_files:
            vals = read_ip_log(f+"_lid.log", test="ldc")
        else:
            vals = read_mixed_log(f+"_lid.log", test="ldc")
#        print "len(vals)", len(vals)
#        print vals
        if f in ["MINI", "CR", "MINI_STAB"]:
            scaled_meshes = [2*m for m in meshes]
            vals = filter(lambda v: v["order"] == 1 and v["n"] in scaled_meshes, vals)
            uniq_vals = []
            for n in scaled_meshes:
                for v in vals:
                    if v["n"] == n:
                        uniq_vals.append(v)
                        break
        else:
            vals = filter(lambda v: v["order"] == 4 and v["n"] in meshes, vals)
            vals.sort(key=lambda v:v["n"])
            uniq_vals = []
            for n in meshes:
                for v in vals:
                    if v["n"] == n:
                        uniq_vals.append(v)
                        break
        div_dict[f] = [v["div"] for v in uniq_vals]
    bargraph_from_dict(div_dict, meshes, log=True, save="div_4_test.pdf",
                       title="4th Order Lid Driven Cavity Divergence", xlabel="Mesh size(n)",
                       ylabel="Divergence Error ($L^2$)", legend=True, loc='center left')



def bargraph_from_dict(h_dict, group_labels, log=False, colors=None, save=None,
                       title=None, xlabel=None, ylabel=None, legend=False, loc='best'):
    bars = len(h_dict) + 2 # Add one for a space between groups
    groups = len(group_labels)

    print h_dict
    for v in h_dict.values():
        if len(v) != groups:
            raise ValueError("h_dict must contain sequences with the same number as labels")

    plt.cla()
    plts = []
    if log:
        plt.yscale("log")

    if colors == "greyscale":
        colors = [str(i*(.75/len(h_dict))) for i in xrange(len(h_dict))]
    elif colors is None or colors == "spectrum":
        colors = "rygbcmk"

    keys = sorted(h_dict.keys())
    for n, k in enumerate(keys):
        plts.append(plt.bar(range(n+1, bars*groups+1, bars), h_dict[k],
                            color=colors[n % len(colors)]))
    plt.xticks([1 + i*bars + bars/2 for i in xrange(groups)], group_labels)

    if log:
        ymin, ymax = plt.ylim()
        ymin = min(reduce(lambda acc, x: acc+x, h_dict.values(), []))/2.0
        plt.ylim((ymin, ymax))
    if legend:
        plt.legend( [p[0] for p in plts], [k.replace('_', '\_') for k in keys], loc=loc)
    if title:
        plt.title(title)
    if ylabel:
        plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    if save is not None:
        plt.savefig(save)



def test_bargraph():
    values = { "CD": [5.78e-4, 9.34e-5, 2.94e-6, 9.25e-8, 2.9e-9],
           "TH" : [2.93e-3, 1.47e-4, 8.16e-6, 4.98e-7, 3.11e-8],
           "STAB" :  [5.78e-4, 9.34e-5, 2.94e-6, 9.25e-8, 2.9e-9],
           "TH_STAB" : [2.93e-3, 1.47e-4, 8.16e-6, 4.98e-7, 3.11e-8],
           "MINI" : [2.93e-3, 1.47e-4, 8.16e-6, 4.98e-7, 3.11e-8],
           "MINI_STAB" : [2.93e-3, 1.47e-4, 8.16e-6, 4.98e-7, 3.11e-8],
           "BF" :  [5.78e-4, 9.34e-5, 2.94e-6, 9.25e-8, 2.9e-9],
           "IP" :  [5.78e-4, 9.34e-5, 2.94e-6, 9.25e-8, 2.9e-9]
           }

    meshes = [4, 8, 16, 32, 64]

    bargraph_from_dict(values, meshes, log=True, colors="spectrum",
                       title="Test Graph", xlabel="mesh size", ylabel="test")

if __name__ == "__main__":
    create_4th_order_plots()
#    test_bargraph()
