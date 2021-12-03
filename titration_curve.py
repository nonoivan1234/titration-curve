import matplotlib.pyplot as plt
import math

'''
Construct a curve for the titration of 25.00mL if 0.1000 M maleic acid, HOOC-CH=CH-COOH, with 0.1000 M NaOH.
We can write the two dissociation equilibria as 
H2M + H2O = H3O+ + HM-  Ka1 = 1.3 * 10^-2
HM- + H2O = H3O+ + M2-  Ka2 = 5.9 * 10^-7
'''

def solve_quadratic_equation(a, b, c):  # 解二次方程式，輸入ax^2 + bx + c = 0 的 a,b,c
    return ((-b+math.sqrt(b**2 - 4*a*c))/(2*a))    # return the greater one root

# Constant
ka1 = 1.3e-2
ka2 = 5.9e-7
kw = 1e-14
kb1 = kw/ka1
kb2 = kw/ka2

# Initial variaties
vol_array = []
pH_array = []
vol_perDrop = 0.1
c0_H2M = 0.1
c0_NaOH = 0.1
vol0_H2M = 25.00
vol_NaOH = 0
pH = 0


# Initial pH
pH = -math.log(solve_quadratic_equation(1, ka1, -ka1*c0_H2M), 10)
vol_array.append(vol_NaOH)
pH_array.append(pH)

vol_NaOH = vol_NaOH + vol_perDrop

# First Buffer Region
'''
[HM-] = c_NaHM + [H3O+] - [OH-]
[H2M] = c_H2M - [H3O+] - [OH-]
Because the solution is quite acidic, [H3O+]>>[OH-]
By ka1 = [H3O+]*[HM-]/[H2M]
We get [H3O+]^2 + (ka1+c_NaHM)*[H3O+] - ka1*c_H2M = 0
Then we can calculate the pH value
'''
while(vol_NaOH < vol0_H2M - vol_perDrop):
    c_NaHM = vol_NaOH*c0_NaOH/(vol0_H2M+vol_NaOH)  
    c_H2M = (vol0_H2M*c0_H2M - vol_NaOH*c0_NaOH)/(vol0_H2M+vol_NaOH)
    
    pH = -math.log(solve_quadratic_equation(1, ka1+c_NaHM, -ka1*c_H2M), 10)
    
    vol_array.append(vol_NaOH)
    pH_array.append(pH)
    vol_NaOH = vol_NaOH + vol_perDrop

# Just Prior to First Equivalence Point
c_HM = vol_NaOH*c0_NaOH/(vol_NaOH+vol0_H2M)
c_H2M = (vol0_H2M*c0_H2M - vol_NaOH*c0_NaOH)/(vol_NaOH+vol0_H2M)
pH = -math.log(solve_quadratic_equation(1+c_HM/ka1, -c_H2M, -ka2*c_HM),10)
vol_array.append(vol_NaOH)
pH_array.append(pH)
vol_NaOH = vol_NaOH + vol_perDrop

# First Equivalence Point 
c_HM = vol_NaOH*c0_NaOH/(vol_NaOH+vol0_H2M)
pH = -math.log(math.sqrt(ka2*c_HM/(1+c_HM/ka1)), 10)
vol_array.append(vol_NaOH)
pH_array.append(pH)
vol_NaOH = vol_NaOH + vol_perDrop

# Just after the First Equivalence Point
while(vol_NaOH < 2*vol0_H2M - vol_perDrop):
    c_NaHM = (vol0_H2M*c0_H2M - (vol_NaOH - vol0_H2M)*c0_NaOH)/(vol0_H2M+vol_NaOH)  
    c_Na2M = (vol_NaOH - vol0_H2M)*c0_NaOH/(vol0_H2M+vol_NaOH)
    
    pH = -math.log(solve_quadratic_equation((ka1 + c_NaHM), ka1*c_Na2M, -ka1*ka2*c_NaHM), 10)
    
    vol_array.append(vol_NaOH)
    pH_array.append(pH)
    
    vol_NaOH = vol_NaOH + vol_perDrop

# Just Prior to Second Equivalence Point
c_NaHM = (vol0_H2M*c0_H2M - (vol_NaOH - vol0_H2M)*c0_NaOH)/(vol0_H2M+vol_NaOH)  
c_Na2M = (vol_NaOH - vol0_H2M)*c0_NaOH/(vol0_H2M+vol_NaOH)
pH = -math.log(kw/solve_quadratic_equation(1, c_NaHM + kb2, -c_Na2M*kb2),10)
vol_array.append(vol_NaOH)
pH_array.append(pH)

vol_NaOH = vol_NaOH + vol_perDrop

# Second Equivalence Point 
c_Na2M = (vol_NaOH - vol0_H2M)*c0_NaOH/(vol0_H2M+vol_NaOH)
pH = -math.log(kw/(math.sqrt(c_Na2M*kb2)), 10)
vol_array.append(vol_NaOH)
pH_array.append(pH)

vol_NaOH = vol_NaOH + vol_perDrop

# pH Just beyond Second Equivalence Point
c_M = vol0_H2M*c0_H2M/(vol0_H2M+vol_NaOH)
c_OH = (vol_NaOH - vol0_H2M*2)*c0_NaOH/(vol0_H2M+vol_NaOH)
c_HM = solve_quadratic_equation(1, c_OH + kb1, -c_M*kb1)
pH = -math.log(kw/(c_OH+c_HM),10)
vol_array.append(vol_NaOH)
pH_array.append(pH)

vol_NaOH = vol_NaOH + vol_perDrop

# pH beyond the Second Equivalence Point
while(vol_NaOH < vol0_H2M*2 + 10):
    c_OH = (vol_NaOH - vol0_H2M*2)*c0_NaOH/(vol0_H2M+vol_NaOH)
    pH = -math.log(kw/c_OH,10)
    vol_array.append(vol_NaOH)
    pH_array.append(pH) 
    
    vol_NaOH = vol_NaOH + vol_perDrop

# Final plot
p = plt.plot(vol_array, pH_array)
plt.title('Titration curve for 25.00 mL 0.1000 M Maleic Acid And 0.1000 M NaOH')
plt.ylabel('pH')
plt.xlabel('Volume of 0.1000 M NaOH, mL')
plt.grid(True)
plt.show()