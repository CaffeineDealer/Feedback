function Fr = normZ(fr,b)


Fr = (exp(-b*fr) * fr) / (0.5 + (5 * b));