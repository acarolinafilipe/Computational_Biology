
There are three functions: Viterbi, forward and backward.

- VITERBI:
Receives an input DNA sequence X and returns the sequence of states that maximizes probability P(x|pi) for all pi paths.
Ex:
Input:
	P=viterbi('CATGCGGGTTATAAC')

Output:

	P =

     2     1     1     2     2     2     2     2     1     1     1     1     1     1     2 

- FORWARD:
Recieves a DNA sequence as input and returns the total probability of the sequence being generated by the Hidden Markov model, as well as the dynamic programming matrix containing forward probabilities. 

Ex:
Input:
	[p1,f]=forward('CATGCGGGTTATAAC')

Output:

	p1 =

   	9.3865e-10


	f =

  	Columns 1 through 13
	
	     0    0.0233    0.0074    0.0019         0    0.0001    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
   	0.1333    0.0092    0.0022    0.0022    0.0007    0.0002    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
    	0.1000    0.0333    0.0057         0    0.0002         0         0         0    0.0000    0.0000    0.0000    0.0000    0.0000

  	Columns 14 through 15

    	0.0000         0
    	0.0000    0.0000
    	0.0000    0.0000 

- BACKWARD:
Recieves a DNA sequence as input and returns the total probability of the sequence being generated by the Hidden Markov model, as well as the dynamic programming matrix containing backward probabilities. 

Ex:
Input:

	[p2,b]=backward('CATGCGGGTTATAAC')

Output:
    	p2 =

   	9.3865e-10


	b =

  	Columns 1 through 13

    	0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0002    0.0009    0.0031    0.0140    0.0494
    	0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0002    0.0009    0.0034    0.0150    0.0548
    	0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0001    0.0003    0.0013    0.0046    0.0209    0.0729

  	Columns 14 through 15

    	0.1600    1.0000
    	0.2750    1.0000
    	0.2500    1.0000

- POSTERIOR PROBABILITIES:
Forward and backward functions can also be used simultaneously to obtain posterior probabilities P(pi_i=k|x). 

Ex:
To obtain P(pi_2=1| x), with x being the input sequence CATGCGGGTTATAAC:

Input:
	[p1,f]=forward('CATGCGGGTTATAAC');
	[p2,b]=backward('CATGCGGGTTATAAC');

	(f(1,2)*b(1,2))/p1

Output:

	ans =
	
	0.404872262414548






