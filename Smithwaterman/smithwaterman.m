function [score, all_align] = smithwaterman(s1,s2,d)
%s1 and s2 are the aminoacid sequences we want aligned;
%d is the module of gap penalty;
%score returns the score of the optimal local alignment(s);
%all_align is a 2xn matrix in which n is the number of optimal alignments

all_align=[]; 

% standardization of the format of the sequences to capital letters
s1 = upper(s1);
s2 = upper(s2);

size1 = length(s1)+1;
size2 = length(s2) +1;

%blosum scores matrix using bioinformatics toolbox
[blosum50, matrix_info] = blosum(50);
blosum50=blosum50(1:20,1:20);

%aminoacid order
amino_order = matrix_info.Order;
amino_order = amino_order(1:20);

%initiliazing scoring matrix
matrixsw = zeros(size2,size1);

%initializing arrow matrix.
%Each position (i,j) of this matrix contains a char sequence indicating the direction(s) of
%the arrow(s) pointing out of element (i,j) in scoring matrix. The char indicates
%an existing arrow: d stands for diagnoal, v for vertical, h for horizontal.
matrixarrows = cell(size2,size1);

%Filling scoring and arrows matrices
for i = 2:size2
    score_i = find(amino_order==s2(i-1));
    
    for j = 2:size1
        score_j = find(amino_order==s1(j-1));
        
        F=[0,matrixsw(i-1,j-1)+blosum50(score_i,score_j),matrixsw(i-1,j)-d,matrixsw(i,j-1)-d];
        
        f_max = find(F==max(F));
        
        %adds the corresponding arrow in char format
        %(diagonal,vertical,horizontal)
        
        if find(f_max==2)
            matrixarrows{i,j}=[matrixarrows{i,j},'d'];
        end
        
        if find(f_max==3)
            matrixarrows{i,j}=[matrixarrows{i,j},'v'];
        end
        
        if find(f_max==4)
            matrixarrows{i,j}=[matrixarrows{i,j},'h'];
        end

        matrixsw(i,j)=F(f_max(1));
        
    end
    
end

%Maximum Score on the Matrix
score=max(matrixsw(:));

%Maximum Score indexes 
[max_i,max_j]=find(matrixsw==score);

if score~=0
    %computing all optimal local alignments
    for i=1:length(max_i)

        c_i = max_i(i); % i index for the ith local alignment 
        c_j = max_j(i); % j index for the ith local alignment 

        seq_al1 = ''; % aligned s1 sequence
        seq_al2 = ''; % aligned s2 sequence 

        all_align=[all_align, find_path(matrixsw,matrixarrows,s1,s2,seq_al1,seq_al2,c_i,c_j)];

    end
end

end