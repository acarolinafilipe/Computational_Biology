%Implementation of viterbi algorithm according to given HMM

%Receives sequence X and, according to HMM model described (here,
%parameters are fixed), computes a path pi* that maximizes P(X,pi) over all 
%possible paths pi

function [pi_star] = viterbi(x) 
%Input:
    %X - character vector of nucleotides, representing a DNA sequence;

%Output:
    %pi_star - array with the same length as X, representing the most probable
%sequence of states for the inputed sequence;

if ~ischar(x)
     error('Please insert valid sequence (it should be a char). Use '' ');
end

%Size of DNA sequence
N = length(x);

%Number of states
k=3;

%Transition matrix
a_11 = 0.6; a_12 = 0.4; a_13 = 0;
a_21 = 0.25; a_22 = 0.5; a_23 = 0.25;
a_31 = 0.25; a_32 = 0.25; a_33 = 0.5;

A = [a_11 a_12 a_13; a_21 a_22 a_23; a_31 a_32 a_33];

a=log(A);
   
%Emission matrix
order='ATGC';

e_1A = 0.4; e_1T = 0.3; e_1G = 0.3; e_1C = 0;
e_2A = 0.1; e_2T = 0.1; e_2G = 0.4; e_2C = 0.4;
e_3A = 0.4; e_3T = 0.3; e_3G = 0; e_3C = 0.3;

E = [e_1A e_1T e_1G e_1C; e_2A e_2T e_2G e_2C; e_3A e_3T e_3G e_3C];

e=log(E);

%Initializing viterbi and ptr matrix
%Both are kxN matrixes
V=zeros(k,N);
ptr=zeros(k,N);

%Since initial probabilities for all the three states are equal, we can
%simply calculate the first column of V with the emission matrix times the
%probability ostates (1/number of states)
V(:,1)=e(:, strfind(order,x(1)));

%Construction of dynamic programming matrix V according to viterbi
%algorithm's formula. Arrows matrix saves the k-1 state that precedes the
%current k state in the path.

for i=2:N
   current_n = strfind(order,x(i));
   for l=1:k
        [m,j] = max(V(:,i-1)+a(:,l));
        V(l,i)= e(l, current_n) + m;
        ptr(l,i) = j;
   end
    
end


%The final state of the path with maximum P(x|pi) is equal to the maximum value of
%the last column of matrix V because transitions from all the states to an end state are equal
[~ , mstate] = max(V(:,N));

%Initializing optimum path
pi_star = zeros(1,N);
pi_star(N)=mstate;

%Traceback of optimum path

i=N;

while(i>1)
    
    pi_star(i-1) = ptr(pi_star(i),i);
     
    i = i-1;
    
end

disp(V)
disp(ptr)
end