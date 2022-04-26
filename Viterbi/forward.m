%Implementation of forward algorithm according to same HMM as viterbi function.

function [P,f] = forward(x) 

%Implementation is similar to the latter, differing in that the
%dynamic programming matrix saves the forward probabilities instead and the path isn't saved. 
%P is the total probability of sequence x being generated by the HMM
%f is the dynamic programming matrix containing forward probabilities

%size of DNA sequence
N = length(x);

%number of states
Q=3;

%transition matrix
a_11 = 0.6; a_12 = 0.4; a_13 = 0;
a_21 = 0.25; a_22 = 0.5; a_23 = 0.25;
a_31 = 0.25; a_32 = 0.25; a_33 = 0.5;

a = [a_11 a_12 a_13; a_21 a_22 a_23; a_31 a_32 a_33];


%emission matrix
%columns correspond to nucleotides while rows correspond to hidden states

order='ATGC';

e_1A = 0.4; e_1T = 0.3; e_1G = 0.3; e_1C = 0;
e_2A = 0.1; e_2T = 0.1; e_2G = 0.4; e_2C = 0.4;
e_3A = 0.4; e_3T = 0.3; e_3G = 0; e_3C = 0.3;

e = [e_1A e_1T e_1G e_1C; e_2A e_2T e_2G e_2C; e_3A e_3T e_3G e_3C];


%initializing f matrix
f=zeros(Q,N);


f(:,1)=e(:, strfind(order,x(1))).*(1/Q);


%Recursive loop computing dynamic programming matrix according to forward
%algorithm's formula

for i=2:N
   current_n = strfind(order,x(i));
   for k=1:Q
        f(k,i)= e(k, current_n)*(f(:,i-1)'*a(:,k));
   end
   
end

%Total probability of x sequence being generated by the model is equal to 
%the sum of all entries of the last column
P=sum(f(:,N));


end