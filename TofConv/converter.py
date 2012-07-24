def input2energy(inputval, option,Ltot=1, theta=1):
    e2lam = 81.787
    e2nu = 4.139
    e2v = 0.0000052276
    e2k = 2.717
    e2t = 0.086165
    ed2cm = 0.123975

    pi = 3.14159265358979323846264
    iv2 = inputval ** 2 


    if option == 'Wavelength':
        Energy = e2lam / iv2

    elif option == 'Energy':
        Energy = inputval
	
	
    elif option == 'Nu':
        Energy = e2nu * inputval

    elif option == 'Velocity':
        Energy = e2v *iv2

    elif option == 'Momentum':
        Energy = e2k*iv2

    elif option == 'Temperature':
        Energy = e2t *inputval

    elif option == 'Invcm':
        Energy = e2cm * inputval

    elif option == 'Q':
        k = inputval * 0.5 / Sin (ThetaRad)
        Energy = e2K * k * k

    elif option == 'Dspacing':
        lam = 2 * inputval * Sin (ThetaRad)
        Energy = e21am / (lam * lam)

    elif  option == 'Tof':
        Energy = 1000000 * Ltot
        Energy = e2v * Energy *Energy / iv2

    return Energy


##---------covert to output units--------------------------------------------------


def energy2output(Energy, option,Ltot=1,theta=1):
    e2lam = 81.787
    e2nu = 4.139
    e2v = 0.0000052276
    e2k = 2.0717
    e2t = 0.086165
    e2cm = 0.123975

    pi = 3.14159265358979323846264
    iv2 = Energy ** 2 

    if option == 'Wavelength':
        OutputVal =  (e2lam/ Energy)**0.5

    elif option == 'Nu':
        OutputVal = Energy / e2nu

    elif option == 'Velocity':
        OutputVal = (Energy / e2v)**0.5

    elif option == 'Momentum':
        OutputVal = (Energy / e2k)**0.5

    elif option == 'Temperature':
        OutputVal = Energy / e2t

    elif option == 'Vcm':
        OutputVal = Energy / e2cm

    elif option == 'Q':
        k = Sqr(Energy / e2k)
        OutputVal = 2 * k * Sin(ThetaRad)

    elif option == 'Dspacing':
        lam = Sqr(e2lam / Energy)
        OutputVal = lam * 0.5 / Sin(ThetaRad)

    elif option == 'Tof':
        OutputVal = Ltot * 1000 * Sqr(e2v * 1000000 / Energy)
    
    elif option == 'Energy':
        OutputVal = Energy
      
    return OutputVal

out = input2energy(10,'Energy',Ltot=10)  
out = energy2output(out, 'Tof')
print(out)

