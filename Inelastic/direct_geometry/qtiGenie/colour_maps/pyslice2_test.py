execfile('PySlice2.py')
ww=data2D()
w1=mtd['w1']
ww.rebinproj(w1,'0,.1,12')
ww.display(10)
ww.CutAlongE(1,5,-10,1,80)
ww.CutAlongE(1,10,-10,1,80,over=True)
ww.CutAlongQ(-10,10,0,.2,12)
ww.CutAlongQ(15,30,0,.2,12,over=True)
