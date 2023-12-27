*double balanced mixer

Vdd 101 0 dc 3
Rload1 101 102 300
Rload2 101 103 300

* mosfets
M1 107 104 102 n 30e-6 .25e-6 2
M2 107 106 103 n 30e-6 .25e-6 2
M3 108 106 102 n 30e-6 .25e-6 2
M4 108 104 103 n 30e-6 .25e-6 2

M5 114 110 107 n 30e-6 .25e-6 2
M6 115 111 108 n 30e-6 .25e-6 2

*source degeneration
*Lde1 114 129 1e-9
Rloss1 114 109 1.2
*Lde2 115 139 1e-9
Rloss2 115 109 1.2

* LC tank
*Lde3 109 149 3e-9
Rloss3 109 0 3.6
*Cde  109 0 9.2e-12

*input
*Vlo+ 154 0 SIN 1 0.6 900e6 0
Vlo+ 154 0 DC 3
Rlo1 154 104 50
*Vlo- 164 0 SIN 1 0.6 900e6 180
Vlo- 164 0 DC -3
Rlo2 164 106 50

*Vrf1+ 112 212 SIN 0.6 0.01 800e6 180
*Vrf2+ 212 0 SIN 0  0.01 600e6 180
Vrf2+ 212 0 dc 3
*Vrf1- 113 213 SIN 0.6  0.01 800e6 0
*Vrf2- 213 0 SIN 0 0.01 600e6 0
Vrf2- 213 0 dc -3
Rs1 112 110 25
Rs2 113 111 25

* level 1 models
.MODEL 1 VT -0.58281 MU 1.224952e-2 COX 6.058e-3 LAMBDA 0.05 CJ0 4.0e-14
.MODEL 2 VT 0.386 MU 3.0238e-2 COX 6.058e-3 LAMBDA 0.05 CJ0 4.0e-14

.hb 100e6 50
.plotnv  102
.plotnv  103
*.plotnv  112
*.plotnv  113
*.plotnv  154
*.plotnv  164
.end