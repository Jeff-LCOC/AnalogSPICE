VDD 2 0 DC 3
Vin 1 0 DC 0

M1 3 1 2 p 30e-6 0.35e-6 1
M2 3 1 0 n 10e-6 0.35e-6 2

.MODEL 1 VT -0.75 MU 5e-2 COX 0.3e-4 LAMBDA 0.05 CJ0 4.0e-14
.MODEL 2 VT 0.83 MU 1.5e-1 COX 0.3e-4 LAMBDA 0.05 CJ0 4.0e-14

.PLOTNV 1
.PLOTNV 3
