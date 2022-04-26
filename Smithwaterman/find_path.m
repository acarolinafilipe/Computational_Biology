function local_al = find_path(matrixsw,matrixarrows,s1,s2,seq_al1,seq_al2,c_i,c_j)
%computes the optimal alignments starting in position given by indexes c_i and c_j;
%matrixsw and matrixarrows are the scoring and arrow matrices, respectively, computed with swmithwaterman function;
%s1 and s2 are the given sequences to align;
%seq_al1 and seq_al2 are the optimal alignments incomplete;


local_al=[]; % Empty array that will store the optimal alignments 

%Traceback of optimal alignment. Stops when first zero in scoring matrix is
%encountered.
while matrixsw(c_i,c_j)~=0
    
    current_arrow = matrixarrows{c_i,c_j};
    
    
    %recursive loop that runs in case of multiple arrows (bifurcation) 
    while length(current_arrow)>1
        new_arrow = current_arrow(1);
        current_arrow=current_arrow(2:end);
        new_ma = matrixarrows;
        new_ma(c_i,c_j)={new_arrow};
        local_al=[local_al, path_find(matrixsw,new_ma,s1,s2,seq_al1,seq_al2,c_i,c_j)];
    end
    
    %extension of alignment based on direction of current arrow
    if current_arrow(1,1)=='d'
        c_i=c_i-1;
        c_j=c_j-1;
        seq_al1=[s1(c_j),seq_al1];
        seq_al2=[s2(c_i),seq_al2];
    
    
    elseif current_arrow(1,1)=='v'
        c_i=c_i-1;
        seq_al2=[s2(c_i),seq_al2];
        
        seq_al1=['-',seq_al1];

        
    elseif current_arrow(1,1)=='h'
        c_j=c_j-1;
        seq_al1=[s1(c_j),seq_al1];
        
        seq_al2=['-',seq_al2]; 
    end
    
    
end

local_al=[local_al, [convertCharsToStrings(seq_al1);convertCharsToStrings(seq_al2)]]; %storing alignment 

end
