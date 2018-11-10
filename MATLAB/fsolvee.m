clear all; close all; clc;

s0 = [0.5;0.5;0.5]
s = fsolve(@nlsistema, s0)