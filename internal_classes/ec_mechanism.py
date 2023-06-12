import numpy as np
from matplotlib import pyplot as plt
from math import exp, log
import re
import os
from helpers import *


class ec_reaction:
    def __init__(
        self,
        edft,
        zpe,
        reactants,
        labels,
        elchem_steps,
        electrode,
        eq_pot,
        tds=None,
        dg=None,
        symfac=1,
        refel="RHE",
        dirlabels=None,
    ):
        
        # placeholder reminder that multiple checks are necessary
        if len({len(edft), len(zpe), len(labels)-1, len(elchem_steps)}) != 1:
            raise ValueError("Mismatch in number of supplied elements")
        
        self.edft = edft
        self.zpe = zpe
        self.tds = [0] * len(self.edft) if tds is None else tds
        self.reactants = reactants
        self.labels = labels
        self.elchem_steps = elchem_steps
        self.electrode = electrode
        self.eq_pot = eq_pot
        
        self.symfac = symfac
        self.refel = refel
        
        if dg is None:
            self.dg = self._calc_dg_zero(self.edft, self.zpe, self.tds, self.reactants, self.symfac)
        else:
            self.dg = dg
            
    
    @staticmethod
    def read_edft(outcar_path):
        for line in read_reverse_order(outcar_path):
            if re.search("sigma", line):
                edft_out = float(line.split()[-1])
                break
        return edft_out

    @staticmethod
    def read_zpe_tds(outcar_path, calc_tds=True):
        kbt = 25.7 #meV
        with open(outcar_path) as outcar:
            zpe_energy = 0.
            ts_energy = 0.
            for line in outcar:
                if re.search("f\s*=", line):
                    zpe_energy += float(line.split()[-2])
                    #b N_a = 1 since we are calculating TS in meV, not kCal/mol
                    ts_energy += kbt * ( (zpe_energy/kbt) / (exp(zpe_energy/kbt)-1) - log(1-exp(-zpe_energy/kbt)) )
            zpe_energy = zpe_energy / 2000 # meV to eV, harmonic approx
            ts_energy = ts_energy / 1000 
            if not calc_tds:
                return zpe_energy
            return (zpe_energy, ts_energy)
            

    @classmethod
    def auto_read(cls, wd, dirlabels, calc_tds=True):
        # automatic read, not fool-proof, assumes /<wd>/<dirlabel>/zpe directory
        edft = []
        zpe = []
        tds = []
        for label in dirlabels:
            edft.append(cls.read_edft(os.path.join(wd, label, "OUTCAR")))
            zpe_energy, ts_energy = cls.read_zpe_tds(os.path.join(wd, label, "zpe", "OUTCAR"), calc_tds)
            zpe.append(zpe_energy)
            tds.append(ts_energy)
        return cls(edft = edft, zpe = zpe, tds = tds)

    
    @classmethod
    def from_dg(cls, dg):
        # probably not necessary
        pass
    
    def _sum_reactants(self, reactants, idx, energy_type):
        out = 0
        for key in reactants.keys():
            out += reactants[key]["reac_part"][idx] * reactants[key]["energies"][energy_type]
        return out
            
    
    def _calc_dg_zero(self, edft, zpe, tds, reactants, symfac):
        dg = []
        for i in range(1, len(edft)):
            de = edft[i] - edft[i-1] + symfac * self._sum_reactants(reactants, i-1, 0)
            dzpe = zpe[i] - zpe[i-1] + symfac * self._sum_reactants(reactants, i-1, 1)
            dtds = tds[i] - tds[i-1] + symfac * self._sum_reactants(reactants, i-1, 2)
            dg.append((de + dzpe - dtds) / symfac)
                                
        dg.append(self.eq_pot*sum(self.elchem_steps) - sum(dg))
        
        return dg
    
    def calc_gmax(self, op):
        dg = np.array(self.dg)

        dg = dg + self.electrode * (self.eq_pot + op) * np.array(self.elchem_steps)
        
        gmax = []
        for i in range(len(dg)):
            gmax.append(np.cumsum(dg[i:]))

        gmax_arr = np.concatenate(gmax)

        position_matrix = []
        for i in range(len(dg)):
            for j in range(len(dg)-i):
                position_matrix.append([i, j+i+1])

        gmax = max(gmax_arr)
        
        gmax_position = position_matrix[np.argmax(gmax_arr)]
        
        return (gmax, gmax_position)
    
    
    def calc_eta_td(self):
        dg = np.array(self.dg)
        eta_td = max(dg) + self.electrode * (self.eq_pot)
        eta_td_position = np.argmax(dg) + 1
        return (eta_td, eta_td_position)
    
    def _g_plot(self, dg_plot, u, print_labels=False, ax=None, **kwargs):

        ax = plt.gca() if ax is None else ax

        right_side = ax.spines["right"]
        top_side = ax.spines["top"]
        right_side.set_visible(False)
        top_side.set_visible(False)

        ax.set_xticks([])
        ax.set_xlim(0, len(self.labels)+1.5)
        ax.set_ylabel('Free Energy, eV')
        ax.set_xlabel('Reaction Coordinate')

        ax.step(np.array(range(len(dg_plot))), dg_plot, where='post', **kwargs)

        ax.text(
            0.95, -0.01,
            rf'$U_\mathrm{{{self.refel}}}=${u:.2f} V',
            verticalalignment='top',
            horizontalalignment='right',
            transform=ax.transAxes,
            color='k', fontsize=12)

        if print_labels:
            for i in range(len(self.labels)):
                ax.text(
                    i + 0.05,
                    dg_plot[i],
                    self.labels[i],
                    fontdict=None,
                    va='bottom',
                    rotation = 0,
                    size = 12)
                
    def _g_plot_trapezo(
        self,
        dg_plot,
        u,
        print_labels=False,
        ax=None,
        # color='pink',
        line_width=0.85,
        **kwargs
    ):
        
        #some duplication of the code from _g_plot, perhaps can be improved / streamlined
        line_shift = 1-line_width
        
        ax = plt.gca() if ax is None else ax
        dg_plot = dg_plot[:-1]
        for i in range(len(dg_plot)):
            ax.plot(
                [i+line_shift, i+line_width],
                [dg_plot[i], dg_plot[i]],
                # color=color,
                **kwargs)
            
            if i > 0:
                ax.plot(
                    [i-line_shift, i+line_shift],
                    [dg_plot[i-1], dg_plot[i]],
                    linestyle='dotted',
                    # color=color,
                    **kwargs
                )
                
            if print_labels:
                ax.text(
                    i + line_shift+0.05,
                    dg_plot[i],
                    self.labels[i],
                    fontdict=None,
                    va='bottom',
                    rotation = 0,
                    size = 12
                )
                
        right_side = ax.spines["right"]
        top_side = ax.spines["top"]
        right_side.set_visible(False)
        top_side.set_visible(False)

        ax.set_xticks([])
        ax.set_xlim(0, len(self.labels)+1.5)
        ax.set_ylabel('Free Energy, eV')
        ax.set_xlabel('Reaction Coordinate')
        ax.text(
            0.95, -0.01,
            rf'$U_\mathrm{{{self.refel}}}=${u:.2f} V',
            verticalalignment='top',
            horizontalalignment='right',
            transform=ax.transAxes,
            color='k', fontsize=12
        )
    

                
    def _g_transform_plot(self, dg, u):
        dg = np.array(dg)
        dg = dg - u * np.array(self.elchem_steps)
        dg_plot = np.cumsum(dg)
        dg_plot = np.concatenate(([0], dg_plot, [dg_plot[-1]]))
        return dg_plot
                
    def g_plot(self, dg=None, u=0, ax=None, trapezo=False, **kwargs):
        dg = self.dg if dg is None else dg
        ax = plt.gca() if ax is None else ax

        dg_plot = self._g_transform_plot(dg, u)
        
        if trapezo:
            self._g_plot_trapezo(dg_plot, u, **kwargs)
        else:
            self._g_plot(dg_plot, u, **kwargs)

        
    def g_plot_eq(self, ax=None, trapezo=False, **kwargs):
        ax = plt.gca() if ax is None else ax
        
        dg_plot = self._g_transform_plot(self.dg, u=self.eq_pot)
        
        if trapezo:
            self._g_plot_trapezo(dg_plot, u=self.eq_pot, print_labels=True, ax=ax, **kwargs)
        else:
            self._g_plot(dg_plot, u=self.eq_pot, print_labels=True, ax=ax, **kwargs)
        
        eta_td, eta_td_position = self.calc_eta_td()

        ax.text(
            eta_td_position + 0.27,
            (dg_plot[eta_td_position-1] + dg_plot[eta_td_position])/2,
            rf'$\eta_\mathrm{{TD}} = $ {eta_td:.2f} V',
            fontdict=None,
            va='center',
            rotation = 270,
            size=12
        )

        ax.annotate(
            '',
            xy=(eta_td_position + 0.25, dg_plot[eta_td_position]),
            xytext=(eta_td_position + 0.25, dg_plot[eta_td_position-1]),
            arrowprops={'arrowstyle':'<->', 'shrinkA': 0, 'shrinkB': 0}
        )
        
    
    def g_plot_gmax(self, op=0.3, ax=None, custom_coord=None, trapezo=False, **kwargs):
            
        ax = plt.gca() if ax is None else ax
        
        dg_plot = self._g_transform_plot(self.dg, u=(self.eq_pot + op))
        
        if trapezo:
            self._g_plot_trapezo(dg_plot, u=self.eq_pot, print_labels=True, ax=ax, **kwargs)
        else:
            self._g_plot(dg_plot, u=(self.eq_pot+op), print_labels=True, ax=ax, **kwargs)
        
        gmax, pos = self.calc_gmax(op)

        if gmax >= 0:
            xytext_coord = (3, -1*op) if custom_coord is None else custom_coord

            ax.annotate(
                '',
                xy=(pos[1] + 0.25, dg_plot[pos[1]]),
                xytext=(pos[1] + 0.25, dg_plot[pos[0]]),
                arrowprops={'arrowstyle':'<->', 'shrinkA': 0, 'shrinkB': 0}
            )
            
            ax.annotate(
                rf'$G_\mathrm{{max}}$({op:.2f} V) =  {gmax:.2f} eV',
                xy = ( pos[1] + 0.25, (dg_plot[pos[1]] + dg_plot[pos[0]])/2 ),
                xytext = xytext_coord, ### this needs to be optimized
                arrowprops={
                    'arrowstyle': '->',
                    'shrinkA': 3,
                    'shrinkB': 3,
                    'connectionstyle':"arc3,rad=0.4",
                    'ls': '--',
                    'color': 'grey',
                    'alpha': 0.5
                }
            )


def construct_ec_mechanism(reactants, labels, elchem_steps, electrode, eq_pot, **kwargs):
    
    class custom_ec_mechanism(ec_reaction):
        def __init__(
            self,
            edft,
            zpe,
            reactants_energies,
            reactants=reactants,
            labels=labels,
            elchem_steps=elchem_steps,
            electrode=electrode,
            eq_pot=eq_pot,
            **kwargs
        ):
            
            self.reactants = reactants
            for key in self.reactants.keys():
                self.reactants[key]["energies"] = reactants_energies[key]
            
            super().__init__(edft, zpe, reactants, labels, elchem_steps, electrode, eq_pot, **kwargs)
            self._check_reactants()
        
        def _check_reactants(self):
            pass
            
    return custom_ec_mechanism

def construct_ec_own(ec_reaction_class, reactants_energies, **kwargs):
    class ec_own(ec_reaction_class):
        def __init__(self, edft, zpe, reactants_energies=reactants_energies, **kwargs):
            super().__init__(edft, zpe, reactants_energies, **kwargs)
            
    return ec_own

oer_mononuc = construct_ec_mechanism(
    # reactants might need to be coded as class attribute instead of instance attribute
    reactants={
        "H2O": {"reac_part": [-1, 0, -1]},
        "H2": {"reac_part": [0.5, 0.5, 0.5]}
    },
    labels=["M", "M-OH", "M-O", "M-OOH", r"M + O$_2(g)$"],
    elchem_steps=[True, True, True, True],
    electrode=-1,
    eq_pot=1.23
)

oer_bifunc1 = construct_ec_mechanism(
    # reactants might need to be coded as class attribute instead of instance attribute
    reactants={
        "H2O": {"reac_part": [-1, 0, -1]},
        "H2": {"reac_part": [0.5, 0.5, 0.5]}
    },
    labels=[
        r"M + *O$_A$",
        r"M-OH + *O$_A$",
        r"M-O + *O$_A$",
        r"M-OO + *OH$_A$",
        r"M + *O$_A$ + O$_2(g)$"],
    elchem_steps=[True, True, True, True],
    electrode=-1,
    eq_pot=1.23
)

oer_bifunc2 = construct_ec_mechanism(
    # reactants might need to be coded as class attribute instead of instance attribute
    reactants={
        "H2O": {"reac_part": [-1, 0, -1, 0]},
        "H2": {"reac_part": [0.5, 0.5, 0, 0.5]}
    },
    labels=[
        r"M + *O$_A$",
        r"M-OH + *O$_A$",
        r"M-O + *O$_A$",
        r"M-OOH + *OH$_A$",
        r"M-OOH + *O$_A$",
        r"M + *O$_A$ + O$_2(g)$"],
    elchem_steps=[True, True, False, True, True],
    electrode=-1,
    eq_pot=1.23
)

oer_binuc = construct_ec_mechanism(
    # reactants might need to be coded as class attribute instead of instance attribute
    reactants={
        "H2O": {"reac_part": [-1, -1, 0, 0]},
        "H2": {"reac_part": [0.5, 0.5, 0.5, 0.5]}
    },
    labels=[
        r"M + M",
        r"M-OH + M",
        r"M-OH + M-OH",
        r"M-O + M-OH",
        r"M-O + M-O",
        r"M + M + O$_2(g)$"],
    elchem_steps=[True, True, True, True, False],
    electrode=-1,
    eq_pot=1.23
)