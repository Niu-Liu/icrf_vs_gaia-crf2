#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: fig_copy.py
"""
Created on Mon May 13 00:37:57 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os

# -----------------------------  FUNCTIONS -----------------------------
# figfiles = [
#     "median-error_bgt.eps",
#     "error-vs-nobs_bgt.eps",
#     "error-vs-nobs-sou-nor_bgt.eps",
#     "pos-err-vs-dec_bgt.eps",
#     "pos-diff-vs-dec_bgt.eps",
#     "gaia-vlbi-separation_bgt.eps",
#     "rotation_bgt.eps",
#     "glide_bgt.eps",
#     "quadrupole_bgt.eps",
#     "nor-sep-dist_bgt.eps",
#     "post_sep_icrf2_vs_icrf3.eps",
#     "nb-sou-on-g.eps",
#     "rotation-on-g.eps",
#     "glide-on-g.eps",
#     "EM20-on-g.eps"]
figfiles = [
    "fig01.eps",
    "fig02.eps",
    "fig03.eps",
    "fig04.eps",
    "fig05.eps",
    "fig06.eps",
    "fig07.eps",
    "fig08.eps",
    "fig09.eps",
    "fig10.eps",
    "fig11.eps",
    "fig12.eps",
    "fig13.eps",
    "fig14.eps"]

dir1 = "../plots/"
# dir2 = "/Users/Neo/github/Niu-LIU/icrf_gaia-crf2/notes/revised/figs/"
dir2 = "../notes/accepted/"

print("Input directory: ", dir1)
print("Output directory: ", dir2)

# Remove old figures
# os.system("rm {}*".format(dir2))

for fig in figfiles:
    # os.system("cp {}{} {}".format(dir1, fig, dir2))
    os.system("mv {}{} {}".format(dir1, fig, dir2))
    print("Copy {}: done!".format(fig))

# --------------------------------- END --------------------------------
