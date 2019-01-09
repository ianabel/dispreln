#!/usr/bin/python

from sys import argv
import netCDF4
import numpy
from scipy import stats

filename = argv[1]

GS2Output = netCDF4.Dataset(filename,'r', format="NETCDF4")

n_kx = len(GS2Output.dimensions["kx"])
n_ky = len(GS2Output.dimensions["ky"])

t_cut = 50
for i_t in range(len(GS2Output["t"])):
    if GS2Output["t"][i_t] > t_cut:
        it_cut = i_t
        break


t_vals = GS2Output["t"][it_cut:]
print "Averaging over the range t= {} to {}".format(t_vals[0],t_vals[-1])

for i_kx in range(n_kx):
    akx = GS2Output["kx"][i_kx]
    for i_ky in range(n_ky):
        aky = GS2Output["ky"][i_ky]
        log_phi_vals = numpy.log(GS2Output["phi2_by_mode"][it_cut:,i_ky,i_kx])/2.0;
        gamma, blah, r_value, blog, err = stats.linregress( t_vals, log_phi_vals );
        if r_value > 0.99:
            print "Mode at k_y = {}, k_x = {} has growth rate {}".format(aky,akx,gamma)
        else:
            print "Mode at k_y = {}, k_x = {} is not displaying exponential time behaviour on the range t = {} to {}".format(aky,akx,t_vals[0],t_vals[-1])

