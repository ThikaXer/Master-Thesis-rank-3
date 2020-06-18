#!/usr/bin/env python
# coding: utf-8

import numpy as np
from itertools import combinations, permutations
import sympy as sp
import gurobipy as gp
from gurobipy import GRB
from tqdm import tqdm
from multiprocessing import Pool
from pathlib import Path
import shutil



# Index-Dictionary for the $\varphi$s

class Beweis_gurobi:
    def __init__(self, points, collinearities, name):
        self.directory = Path(name)
        if self.directory.is_dir():
            shutil.rmtree(self.directory)
        self.directory.mkdir()

        # Initialise and alternate linear and non-linear steps
        self.points = points
        self.collinearities = collinearities
        self.phis = [comb for comb in combinations(points, 3)
                if set(comb) not in collinearities]  # List of bases != 0
        self.phi_idx = {frozenset(phi): i for i, phi in enumerate(self.phis)}  # Dictionary of all indices from corresponding bases
        self.cross_ratios = set()
        for x in points:
            for a, b, c, d in permutations(points - {x}, 4):
                if {x, a, c} not in self.collinearities \
                        and {x, a, d} not in self.collinearities \
                        and {x, b, c} not in self.collinearities \
                        and {x, b, d} not in self.collinearities:
                    self.cross_ratios.add((a, b, c, d, x))
        self.init_GPR()
        self.remaining_cr = len(self.cross_ratios)

        print('Number of cross ratios to check:', len(self.cross_ratios))
        for i in range(1, 11):
            self.linear_step()
            if len(self.cross_ratios) == 0:
                with open(self.directory/'Status.txt', 'w') as f:
                    print(f'Proof is done after linear step with {i} iterations! All cross-ratios are real.', file=f)
                print(f'Proof is done after linear step with {i} iterations! All cross-ratios are real.')
                break
            print('Number of cross ratios left:', len(self.cross_ratios))
            if self.remaining_cr == len(self.cross_ratios):
                with open(self.directory/'Status.txt', 'w') as f:
                    print(f'No Proof. Iteration stationary with {self.remaining_cr} cross-ratios left.', file=f)
                print(f'No Proof. Iteration stationary with {self.remaining_cr} cross-ratios left.')
                break
            else:
                self.remaining_cr = len(self.cross_ratios)
            self.non_linear_step()

        np.savetxt(self.directory/'GPR_matrix.txt', self.GPR, fmt='%i')



    def get(self, a, b, c):
        return self.phi_idx[frozenset([a, b, c])]  # returns index of basis a,b,c

    # methods to build and show a cross ratio vector

    def cross_ratio_vector(self, a, b, c, d, x):
        # builds cross-ratio-vector with 1 for elements with pos exponent, -1 for elements with neg exponent and 0 else
        gpr = np.zeros(len(self.phis), dtype=np.int)
        gpr[self.get(x, a, c)] = 1
        gpr[self.get(x, b, d)] = 1
        gpr[self.get(x, a, d)] = -1
        gpr[self.get(x, b, c)] = -1
        return gpr  # returns a vector

    def vector_to_equation(self, cr):
        # returns cross ratio from cross-ratio vector
        try:
            i1, i2 = np.argwhere(cr == 1).flat
        except Exception as e:
            print(e)
            print(cr)
        links1 = f'φ{tuple(self.phis[i1])}'
        links2 = f'φ{tuple(self.phis[i2])}'

        j1, j2 = np.argwhere(cr == -1).flat
        rechts1 = f'φ{tuple(self.phis[j1])}'
        rechts2 = f'φ{tuple(self.phis[j2])}'

        equation = links1 + links2 + ' = ±' + rechts1 + rechts2
        return equation

    # All possibilities for cross ratios $cr(a, b \mid x, y)_c$, where $a, b, c$ are collinear

    def init_GPR(self):
        self.GPR = set()  # Matrix of all known real-valued cross-ratios
        for coll in self.collinearities:
            others = self.points - coll  # all possibilities for x,y
            for c in coll:  # all possible viewpoints
                a,b = [l for l in coll if l != c]
                for x, y in combinations(others, 2): # all other values for x,y
                    for a_permuted, b_permuted, x_permuted, y_permuted in permutations([a,b,x,y]):
                        if {c, a_permuted, x_permuted} not in self.collinearities \
                                and {c, a_permuted, y_permuted} not in self.collinearities \
                                and {c, b_permuted, x_permuted} not in self.collinearities \
                                and {c, b_permuted, y_permuted} not in self.collinearities:
                            # check all permutations from S_4 in cross ratio
                            gpr_permuted = self.cross_ratio_vector(a_permuted, b_permuted, x_permuted, y_permuted, c)
                            # if permuted cross-ratio already in matrix GPR,
                            # we do not need to ppend it again
                            positive = tuple(gpr_permuted)
                            negative = tuple(-gpr_permuted)
                            if positive not in self.GPR and negative not in self.GPR:
                                self.GPR.add(max(positive, negative))

        self.GPR = np.stack(list(reversed(sorted(self.GPR))), axis=0)  # returns matrix from sorted list
        print(self.GPR.shape)

    # signs are disregarded as we are only interested in the cross ratio being real or not

    # Calculate if vector of cross-ratio $cr(a, b \mid c, d)_x$ is in the integer span by solving integer linear feasibility problem with gurobi

    def solve_lin_span(self, abcdx):
        a, b, c, d, x = abcdx
        vector = self.cross_ratio_vector(a, b, c, d, x)
        try:
            model = gp.Model("Model for vector")
            z = []
            for i in range(len(self.GPR)):
                z.append(model.addVar(vtype=GRB.INTEGER, lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"z_{i}"))
                # z.append(model.addVar(vtype=GRB.INTEGER, lb=-1, ub=1, name=f"z_{i}"))
            for n in range(len(self.GPR[0])):
                model.addConstr(sum(int(self.GPR[m][n]) * z[m] for m in range(len(self.GPR))),
                    GRB.EQUAL, vector[n], f"constraint_{n}")
            # model.setObjective(sum(z[i]*z[i] for i in range(len(self.GPR))), GRB.MINIMIZE)
            model.setObjective(1)
            model.Params.OutputFlag = 0
            model.Params.DualReductions = 0
            model.Params.MIPFocus = 1
            model.optimize()
            if model.Status == 3: #Status = 3 if System is infeasible
                solution = None
            elif model.Status == 2:
                solution = np.array([var.X for var in z])
            else:
                print("error")
            print(model.Status, end='', flush=True)
        except gp.GurobiError as e:
            print('Error code ' + str(e.errno) + ': ' + str(e))
        return solution

    # Linear step: search cross-ratios, which are in span of matrix GPR andwrite them in list new_cross_ratios
    def linear_step(self):
        self.new_cross_ratios = set()
        GPR = list(self.GPR)

        with Pool() as p:
            results = p.map(self.solve_lin_span, self.cross_ratios)

        print()

        GPR = set(tuple(row) for row in self.GPR)
        for cr, solution in zip(self.cross_ratios, results):
            if solution is not None:
                with open(self.directory/'combos_for_new_cross_ratios.txt', 'a+') as f:
                    a,b,c,d,x = cr
                    print(f'cr({a},{b}|{c},{d})_{x}: φ({a}, {c}, {x})φ({b}, {d}, {x}) = ±φ({a}, {d}, {x})φ({b}, {c}, {x})', file = f)
                    for l in range(len(solution)):
                        if solution[l] != 0:
                            for _ in range(abs(int(solution[l]))):
                                print(self.vector_to_equation(self.GPR[l] * np.sign(solution[l])), file=f)
                    print(file=f)


                self.new_cross_ratios.add(cr)
                gpr = self.cross_ratio_vector(*cr)

                positive = tuple(gpr)
                negative = tuple(-gpr)
                if positive not in GPR and negative not in GPR:
                    GPR.add(max(positive, negative))

        self.cross_ratios = self.cross_ratios - self.new_cross_ratios
        self.GPR = np.stack(list(reversed(sorted(GPR))), axis=0)

    # Non-linear step: build for cross-ratios in new_cross_ratios all permutations from $S_4$ and append corresponding vectors to the matrix GPR

    def non_linear_step(self):
        GPR = set(tuple(row) for row in self.GPR)
        # build all permutations of cross-ratios in new_cross_ratios
        for a, b, c, d, x in tqdm(self.new_cross_ratios, desc = 'Non linear step', ncols=120):
            for pi_a, pi_b, pi_c, pi_d in permutations([a, b, c, d]):
                if {x, pi_a, pi_c} not in self.collinearities \
                        and {x, pi_a, pi_d} not in self.collinearities \
                        and {x, pi_b, pi_c} not in self.collinearities \
                        and {x, pi_b, pi_d} not in self.collinearities:
                    gpr_pi = self.cross_ratio_vector(pi_a, pi_b, pi_c, pi_d, x)

                    positive = tuple(gpr_pi)
                    negative = tuple(-gpr_pi)
                    if positive not in GPR and negative not in GPR:
                        GPR.add(max(positive, negative))

        self.GPR = np.stack(list(reversed(sorted(GPR))), axis=0)  # build matrix from sorted list
