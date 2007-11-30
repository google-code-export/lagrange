import sys
from pprint import pprint
import scipy
from scipy.integrate import odeint
from scipy import optimize
array = scipy.array

class GSE2_Parameters:
    labels = ("lam_a", "lam_b", "lam_ab",
              "mu_a", "mu_b", "q_ab", "q_ba")
    def __init__(self, v=None):
        if v:
            assert len(v) == 7
            for i, s in enumerate(self.labels):
                setattr(self, s, v[i])
        else:
            self.lam_a = 0.3
            self.lam_b = 0.4
            self.lam_ab = 0.9
            self.mu_a = 0.4
            self.mu_b = 0.4
            self.q_ab = 0.1
            self.q_ba = 0.1

    def v(self):
        return tuple([ getattr(self, x) for x in self.labels ])

def dDNa_dt(lam_a, mu_a, q_ab, D_Na, D_Nab, E_a):
    v = -(lam_a + mu_a + q_ab)*D_Na \
        + q_ab*D_Nab \
        + 2*lam_a*E_a*D_Na
    return v

def dDNb_dt(lam_b, mu_b, q_ba, D_Nb, D_Nab, E_b):
    v = -(lam_b + mu_b + q_ba)*D_Nb \
        + q_ba*D_Nab \
        + 2*lam_b*E_b*D_Nb
    return v

def dDNab_dt(lam_a, lam_b, lam_ab, mu_a, mu_b,
             D_Na, D_Nb, D_Nab, E_b, E_a, E_ab):
    v = -(lam_a+lam_b+lam_ab+mu_a+mu_b)*D_Nab \
        + lam_a*E_a*D_Nab \
        + lam_a*E_ab*D_Na \
        + lam_b*E_b*D_Nab \
        + lam_b*E_ab*D_Nb \
        + lam_ab*E_a*D_Nb \
        + lam_ab*E_b*D_Na \
        + mu_a*D_Nb \
        + mu_b*D_Na
    return v

def dEa_dt(mu_a, q_ab, lam_a, E_a, E_ab):
    v = mu_a \
        - (mu_a+q_ab+lam_a)*E_a \
        + q_ab*E_ab \
        + lam_a*(E_a**2)
    return v

def dEb_dt(mu_b, q_ba, lam_b, E_b, E_ab):
    v = mu_b \
        - (mu_b+q_ba+lam_b)*E_b \
        + q_ba*E_ab \
        + lam_b*(E_b**2)
    return v

def dEab_dt(lam_a, lam_b, lam_ab, mu_a, mu_b, E_a, E_b, E_ab):
    v = mu_a*E_b \
        + mu_b*E_a \
        - (mu_a+mu_b+lam_a+lam_b+lam_ab)*E_ab \
        + lam_a*E_ab*E_a \
        + lam_b*E_ab*E_b \
        + lam_ab*E_a*E_b
    return v

def dE_dt(y, t, p):
    "p is a GSE2_Parameters object"
    E_a, E_b, E_ab = y
    v = array([
        dEa_dt(p.mu_a, p.q_ab, p.lam_a, E_a, E_ab),
        dEb_dt(p.mu_b, p.q_ba, p.lam_b, E_b, E_ab),
        dEab_dt(p.lam_a, p.lam_b, p.lam_ab, p.mu_a, p.mu_b, E_a, E_b, E_ab)
        ])
    return v

def dDN_dt(y, t, E_a, E_b, E_ab, p):
    "p is a GSE2_Parameters object"
    D_Na, D_Nb, D_Nab = y
    v = array([
        dDNa_dt(p.lam_a, p.mu_a, p.q_ab, D_Na, D_Nab, E_a),
        dDNb_dt(p.lam_b, p.mu_b, p.q_ba, D_Nb, D_Nab, E_b),
        dDNab_dt(p.lam_a, p.lam_b, p.lam_ab, p.mu_a, p.mu_b,
                 D_Na, D_Nb, D_Nab, E_b, E_a, E_ab)
        ])
    return v

def integrate_branch(DN, t, p):
    """
    given probabilities D_Na, D_Nb, and D_Nab at the end of the branch
    (t0), calculate the probabilities at t0+t (rootward)
    p is a GSE2_Parameters object
    """
    E = odeint(dE_dt, (0.0, 0.0, 0.0), (0.0, t), args=(p,))[1]
    D = odeint(dDN_dt, DN, (0.0, t), args=tuple(list(E)+[p]))[1]
    return D

## p = GSE2_Parameters((0.1, 0.1, 0.3, 0.03, 0.03, 0.1, 0.05))
## print integrate_branch((1,0,0), 0.1, p)
## sys.exit()

def evaluate_node(DN, DM, p):
    D_Na, D_Nb, D_Nab = DN
    D_Ma, D_Mb, D_Mab = DM
    v = (
        D_Na*D_Ma*p.lam_a,
        D_Nb*D_Mb*p.lam_b,
        D_Na*D_Mb*p.lam_ab + D_Nb*D_Ma*p.lam_ab \
        + D_Na*D_Mab*p.lam_a + D_Nb*D_Mab*p.lam_b \
        + D_Nab*D_Ma*p.lam_a + D_Nab*D_Mb*p.lam_b
        )
    return v

p = GSE2_Parameters((0.1, 0.1, 0.3, 0.03, 0.03, 0.05, 0.05))
DN = integrate_branch((0,0,1), 1, p)
DM = integrate_branch((0,0,1), 1, p)
print evaluate_node(DN, DM, p)
sys.exit()

def evaluate_node_downpass(node, data, p):
    "p is a GSE2_Parameters object"
    v = []
    for c in node.children():
        if not c.istip:
            v.append(evaluate_node_downpass(c, data, p))
        else:
            DN = integrate_branch(data[c.label], c.length, p)
            v.append(DN)
    assert len(v) == 2, "polytomy at node %s" % node.label
    DN, DM = v
    rv = evaluate_node(DN, DM, p)
    #print node.label, rv
    return rv

def d_na_dT(lam_a, mu_a, mu_b, n_a, n_ab):
    return lam_a*n_a - mu_a*n_a + mu_b*n_ab

def d_nb_dT(lam_b, mu_a, mu_b, n_b, n_ab):
    return lam_b*n_b - mu_b*n_b + mu_a*n_ab

def d_nab_dT(q_ab, q_ba, n_a, n_b):
    return q_ab*n_a + q_ba*n_b

def dnx_dT(y, t, p):
    "p is a GSE2_Parameters object"
    n_a, n_b, n_ab = y
    v = array([
        d_na_dT(p.lam_a, p.mu_a, p.mu_b, n_a, n_ab),
        d_nb_dT(p.lam_b, p.mu_a, p.mu_b, n_b, n_ab),
        d_nab_dT(p.q_ab, p.q_ba, n_a, n_b)
        ])
    return v

def integrate_nx_forward(n_a, n_b, n_ab, t):
    y = (n_a, n_b, n_ab)
    return odeint(dnx_dT, y, (0.0, t))[1]

#print integrate_nx_forward(1,0,0,10)
#sys.exit()



if __name__ == "__main__":
    import newick, phylo
    s = "((((((((P_hawaiiensis_WaikamoiL1:0.010853,P_mauiensis_Eke:0.010853):0.007964,(P_fauriei2:0.013826,P_hathewayi_1:0.013826):0.004991):0.001986,(P_kaduana_PuuKukuiAS:0.020803,P_mauiensis_PepeAS:0.020803):1e-05):0.003762,P_kaduana_HawaiiLoa:0.024565):0.003398,(P_greenwelliae07:0.01271500,P_greenwelliae907:0.01271500):0.01524800):0.018984,((((P_mariniana_MauiNui:0.02241,P_hawaiiensis_Makaopuhi:0.02241):0.008236,P_mariniana_Oahu:0.030646):0.002893,P_mariniana_Kokee2:0.033539):0.005171,P_wawraeDL7428:0.03871):0.008237):0.008255,(P_grandiflora_Kal2:0.027864,P_hobdyi_Kuia:0.027864):0.027338):0.003229,((P_hexandra_K1:0.026568,P_hexandra_M:0.026568):0.005204,P_hexandra_Oahu:0.031771):0.026659);"
    tree = newick.parse(s)
    phylo.polarize(tree)
    for i, n in enumerate(tree.descendants()):
        if n.parent:
            n.length *= 100
        if not n.istip:
            n.label = "N%s" % i

    data = {
        "P_hawaiiensis_WaikamoiL1":  (1,0,0),
        "P_mauiensis_Eke":           (0,1,0),
        "P_fauriei2":                (0,1,0),
        "P_mariniana_Kokee2":        (1,0,0),
        "P_mariniana_Oahu":          (1,0,0),
        "P_mariniana_MauiNui":       (0,1,0),
        "P_hawaiiensis_Makaopuhi":   (1,0,0),
        "P_wawraeDL7428":            (1,0,0),
        "P_kaduana_PuuKukuiAS":      (0,1,0),
        "P_mauiensis_PepeAS":        (0,0,1),
        "P_hathewayi_1":             (0,1,0),
        "P_kaduana_HawaiiLoa":       (0,1,0),
        "P_greenwelliae07":          (1,0,0),
        "P_greenwelliae907":         (1,0,0),
        "P_grandiflora_Kal2":        (1,0,0),
        "P_hobdyi_Kuia":             (1,0,0),
        "P_hexandra_K1":             (1,0,0),
        "P_hexandra_M":              (1,0,0),
        "P_hexandra_Oahu":           (1,0,0),
        }

    tup2lab = {(1,0,0): "A", (0,1,0): "B", (0,0,1): "AB"}
    d = dict([ (x.label, tup2lab[data[x.label]]) for x in tree.leaves() ])

    def f(params):
        lam_w, lam_b, mu, q = params
        for x in params:
            if (x <= 0) or x >= 10:
                return 10e10
        if (lam_w + lam_b) > 10:
            return 10e10
        p = GSE2_Parameters((1.0, 1.0, lam_b, mu, mu, q, q))
        v = evaluate_node_downpass(tree, data, p)
        rv = -sum(scipy.log(v))
        if (rv < 0) or rv == scipy.nan:
            return 10e10
        return rv

    opt_p, neglnL = optimize.fmin_powell(f, [0.5, 0.5, 0.1, 0.1])

    print tree.render_ascii(scaled=False, minwidth=80, data=d)
