count 300
cset T_head=150
changetitle Vycor +H2O E=1.5KV Cooling 50meV 150Hz 4x4.5
count 200
count 200
count 100
changetitle Vycor +H2O E=1.5KV T=200K 50meV 150Hz 4x4.5
count 800
cset T_head=300
changetitle Vycor +H2O E=1.5KV Heating 50meV 150Hz 4x4.5
count 200
count 200
count 150
changetitle Vycor +H2O E=1.5KV T=300K 500meV 150Hz 4x4.5
count 800
scanTemp 290 -20 200 200 Vycor+H20 E=0V Ei=150mev
setWhiteBeam 17881
reduction_ei 300
defineRoi 1 10 20 200
countStats 4 10