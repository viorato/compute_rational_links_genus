#!/usr/bin/env python
# Minimal genus surfaces for rational links
# Jesus Rodriguez-Viorato,  jesusr@cimat.mx,
# Enrique Ramirez-Losada,  kikis@cimat.mx,
# Fabiola Manjarrez Gutierrez, fabiola.manjarrez@im.unam.mx,
# and 
# Mario Eudave Munoz, eduave@matem.unam.mx
# May, 2019

# This program will compute the mininal genus of a surface in the exterior of a 2-bridge link (K_1 \cup K_2) , 
# with the condition that it has only one longitudinal component on K_1. 
# 
# The present code is based on the code written by Jim Hoste and Patrick D. Shanahan, so many of their function
# are being used.

import sys, re
from fractions import Fraction
from continued_fractions import CFraction
 
out = sys.stdout.write
get_input = sys.stdin.readline
 
################################################## 
def gcd(a,b):
    """ A simple routine to compute gcd, needed to test input in UI
    """
    while b:
        a,b = b, a % b
    return a
 
##################################################
def get_even_contiuned_fraction(a,b):
    """ Find the continued fraction [0, 2a_1, 2a_2, ..., 2a_k] of the rational a/b
    """
    if  b<0:
        b = (-1)*b
        a = (-1)*a
    if b == 1:
        return [a]   
    elif b > a and a > (-1)*b :
        return [0]+get_even_contiuned_fraction(b,a)
    else:
        q, r = divmod(a,b)
        if q%2 ==0:
            return [q]+get_even_contiuned_fraction(b,r)
        else:
            return [q+1]+get_even_contiuned_fraction((-1)*b,b-r)

##################################################         
def ask_question(Question, Answers):
    while 1:
        out(Question)
        inp = get_input()[0]    
        if re.search(inp, Answers):
            break
        out("Sorry, input not understood.\n")
    return inp

##################################################
def sgn(n):
    """ Returns the sign of the number n. If n=0 resturns 0.
    """
    if n > 0:
        return 1
    elif n < 0:
        return -1
    else:
        return 0

##################################################
def expanded_cf(cf):
    """ Computes the expanded form of the continue fraction [a_0, 2a_1, 2a_3, ..., 2a_n]
        by replacing each 2a_i with i odd with a block of the form [2,0,2, .., 0, 2] 
        if a_i > 0 or [-2, 0, -2, ... , 0, -2] if a_i < 0
    """
    count = 0
    expanded = [cf[0]]
    for i in cf[1:]:
        if count %2 == 0:
            expanded += [sgn(i)*2]+[0,sgn(i)*2]*(abs(int(i/2))-1)
        else:
            expanded.append(i)
        count +=1
    return expanded

################################################## 
def reducePath(path):
    """ This funtion is no longer in use. 
        The computed paths are reduced.
    """
    l = len(path)
    if l <= 2:
      return path
    reduce_path = [path[0]]
 
    i=0
    while i < l-2:
        if path[i] != path[i+2]:
          reduce_path+=[path[i+1], path[i+2]]
        i+=2   
 
    if l%2 == 0:
      reduce_path+=[path[l-1]]
 
    return reduce_path

################################################# 
def compute_path0(cf):
    """ Transform a continued fraction cf = [a_0, a_1, a_2, ..., a_{2k+1}] into a set of fractions
        [r_0, r_1, ..., r_k] where r_i = [a_0, ..., a_{2i+1}]
    """
    ecf = expanded_cf(cf)
    path0 = []
    for i in range(1,len(ecf),2):
        fr = CFraction(ecf[0:(i+1)])
        path0.append(fr.fraction())
    return path0

################################################# 
def compute_paths(cf):
    """ Transform a continued fraction cf = [a_0, a_1, a_2, ..., a_n] into a set of fractions
        [r_0, r_1, ..., r_n] where r_i = [a_0, ..., a_i]
    """
    ecf = cf
    path1 = []
    for i in range(len(ecf)):
        fr = CFraction(ecf[0:(i+1)])
        path1.append(Fraction(fr.fraction()))
    return path1

################################################# 
def compute_word(cf):
    """ Convert a continued fraction into a string of A's and D's.
    """
    count=0
    cadena = 'A'
    for i in cf[1:]:
        if count%2 == 0:
           if i!=0:
            cadena += 'D'*(int(abs(i)/2))
        else:
            cadena += 'AA'
        count+=1  
    return cadena+'A'


    
################################################# 
def compute_euler(word, mu):
    """ Compute the euler number of the surface given by a word of A's and D's, 
        we assume rho = 1 and mu is a required parameter.
    """
    euler = mu+1
    for char in word:
        if char=='D':
            euler-= (mu -1)
        else:
            euler-=1
    return euler

#################################################
def path2str(path):
    """ Convert a vector of objects into a string representing it.
        Helpfull for getting nice screen outputs of paths.
    """ 
    string = "[ 1/0 , "
    for p in path:
        string += str(p.numerator)+'/'+str(p.denominator) + ' , '
    return string[:-2]+']'


#################################################
def format_surface_info(name, genus, slope1, slope2, boundaries):
    """ This function gives format to display surfaces info.
    """
    if slope2 is not False:
      return 'Genus({0}) = {1}, slope_1({0}) = {2}, slope_2({0}) = {3}, |boundary({0})| = {4}'.format(name, genus, slope1, slope2, boundaries)
    else:
      return 'Genus({0}) = {1}, slope_1({0}) = {2}, slope_2({0}) = empty, |boundary({0})| = {3}'.format(name, genus, slope1, boundaries)

#################################################
def format_more_info(linking, wrapping, paths_info, level=1):
    """ This function gives format to display extra info, such as linking, wrapping, and computed paths.
    """
    tabs = '\t'*level
    text = '{0}LINKING = {1}\n'.format(tabs, linking)
    text+= '{0}WRAPPING = {1}\n'.format(tabs, wrapping)  
  
    for index, path_info in enumerate(paths_info):
        text += tabs + 'SURFACE_{0} PATH IN D{1}: {2}\n'.format(index, path_info['t'], path_info['word'])    
        text += tabs + '\tPATH_{0}: {1}\n'.format(index, path2str(path_info['path']))
   
    #print('\tCONTINUED FRACTIONS:\n\t\tEven: {0} \n\t\tOne+Even: {1}\n'.format( cf_even, cf_1even))
    return text[:-1]

##########################################################
# This is a simple front-end that asks the user for 
# p and q and then computes the boundary slope of the
# p/q 2-bridge link. Output is reformatted2.30-1 to match
# Table 5 in "Computing boundary slopes of 2-bridge links"
# by Hoste and Shanahan
 
if __name__ == "__main__":# main program
    while 1: # get p/q  and compute boundary slopes until user chooses to quit
        while 1: #get p/q until valid input
            while 1: # get p until valid input
                try:
                    alpha=int(input("Enter alpha. "))
                    break
                except ValueError:
                    print("Not a valid integer. Try again.")
            while 1: # get q until valid input
                try:
                    beta=int(input("Enter beta. "))
                    break
                except ValueError:
                    print ("Not a valid integer. Try again.")

            #see if p and q are valid input
            if 0<alpha and alpha<beta and beta%2==0 and gcd(alpha,beta)==1:
                break
            else:
                print("You must enter 0<alpha<beta, beta even, and gcd(alpha,beta)=1")
        while 1:
            try:
                mu = int(input("Enter mu."))
                break
            except ValueError:
                print ("Not a valid integer. Try again.")


        # Compute the continued fractions with only even numbers
        cf_even = get_even_contiuned_fraction(alpha,beta)
        cf_1even = get_even_contiuned_fraction(alpha-beta,beta)
        cf_1even[0]=1
 
        # Computing linking and wrapping using the continued fraction
        linking = 0
        wrapping = 0
        for i in cf_even[1::2]:
            linking+=i
            wrapping += abs(i)
        linking = linking/2
        wrapping = wrapping/2 

        #Number of boundary components
        n = gcd(mu, linking)

        # Actual paths 
        path1 = compute_paths(cf_even)
        path2 = compute_paths(cf_1even)

        #Compute words
        word_even = compute_word(cf_even)
        word_1even = compute_word(cf_1even)

        #Compute Euler number
        euler_even = compute_euler(word_even, mu)
        euler_1even = compute_euler(word_1even, mu)

        #Compute genus
        genus_even = (2-n-euler_even)/2
        genus_1even = (2-n-euler_1even)/2

        print('')  
        print('COMPUTING MINIMAL GENUS SURFACE FOR:\nRATIONAL LINK {0}/{1} WITH rho = 1 and mu = {2}.'.format(alpha, beta, mu))
        print('')

        #Case mu=0 only has surfaces when linking = 0.
        if mu == 0:
            if linking ==0:
                print('THERE IS ONLY ONE SURFACE F AND IS GIVEN BY A D-TYPE EDGE-PATH IN D0')
                genus = wrapping/2               
                print('\t' +  format_surface_info('F', genus, 0, False, 1))
                
                if ask_question('Do you want more info? (y/n)', 'yn') == 'y':
                    path0 = compute_path0(cf_even)
                    print(format_more_info(linking, wrapping, [{'t':0, 'word':'D'*len(path0), 'path':path0}]))

            else:
                print('THERE ARE NO SURFACES WITH mu=0, BECAUSE linking IS NOT ZERO')
            
                   
        # Case mu=0 has two surfaces given by A-Type edge-path.
        if mu == 1:
            print('THERE ARE TWO SURFACES F_1 AND F_2 GIVEN BY A-TYPE EDGE-PATHS IN D1')

            genus1 = (len(cf_even)-2)/2
            genus2 = (len(cf_1even)-2)/2

            print('\t' +  format_surface_info('F_1', genus1, 0, 0, 2))
            print('\t' +  format_surface_info('F_2', genus2, 0, 0, 2))

            if ask_question('Do you want more info? (y/n)', 'yn') == 'y':
                surfaces_info = [
                    {'t':1, 'word':'A'*len(cf_even),'path':path1},
                    {'t':1, 'word':'A'*len(cf_1even),'path':path2}
                ]
                print(format_more_info(linking, wrapping, surfaces_info))

        #Case mu>1 is the trickiest one, it could have AB- and AD-edge-paths.
        if mu>1:

            # Checking if the even path could make a connected surface and if 
            # it is minimal
            check_even = True
            if(len(cf_even) > 3):
                for a in cf_even[1:]:
                    if abs(a) <= 2:
                       check_even = False
                       break
            else:
                check_even = False

            # Checking if the one+even path could make a connected surface and if 
            # it is minimal
            check_1even = True        
            if(len(cf_1even) > 3):
                for a in cf_1even[1:]:
                    if abs(a) <= 2:
                       check_1even = False
                       break
            else:
                check_1even = False

            paths_info = []
            if check_even or check_1even:
                print("THERE ARE SURFACES GIVEN BY AB-TYPE EDGEPATHS IN Dt:")
                
            if check_even:
                genus = (len(cf_even)-4)*(n+1)/4 + 1
                print('\t' + format_surface_info('F_1', genus, 0, 0, n+1))
                paths_info.append({'t':Fraction(1,mu), 'word':'ABBA'*((len(cf_even)-2)/2), 'path': path1})
            if check_1even:
                genus = (len(cf_1even)-4)*(n+1)/4 + 1
                print('\t' + format_surface_info('F_1', genus, 0, 0, n+1))
                paths_info.append({'t':Fraction(1,mu), 'word':'ABBA'*((len(cf_1even)-2)/2), 'path': path2})

            print("THERE ARE TWO SURFACES GIVEN BY AD-TYPE EDGEPATHS IN Dt:")
            print('\t' + format_surface_info('F_1', genus_even, Fraction(-1*linking,mu), Fraction(1*linking*mu,1), n+1))
            print('\t' + format_surface_info('F_2', genus_1even, Fraction(linking,mu), Fraction(linking*mu,1), n+1))

            paths_info.append({'t':Fraction(1,mu), 'word':word_even, 'path': path1})
            paths_info.append({'t':Fraction(1,mu), 'word':word_1even, 'path': path2})

            if ask_question('Do you want more info? (y/n)', 'yn') == 'y':
                print(format_more_info(linking, wrapping, paths_info))


        if ask_question('Do you want to compute another? (y/n)', 'yn') == 'n':
            break
        
