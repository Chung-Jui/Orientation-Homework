function [Q,sigma,R] = FEA(r,F,E,L)
    node_coor = L.*[2 1; 2 0; 1 1; 1 0; 0 1; 0 0];
    element_table = [3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4];
    
    %derive cross-sectional area, the lengths of each element and the cosine and sine of their orientation.
    for i=1:10
        element_table(i,3) = E;
        if i<=6
            element_table(i,4) = pi*r(1)*r(1);
        end
        if i>6
            element_table(i,4) = pi*r(2)*r(2);
        end
        x_diff = (node_coor(element_table(i,2),1)- node_coor(element_table(i,1),1));
        y_diff = (node_coor(element_table(i,2),2)- node_coor(element_table(i,1),2));
        element_table(i,5) = sqrt(x_diff^2 + y_diff^2);
        element_table(i,6) = x_diff/element_table(i,5);
        element_table(i,7) = y_diff/element_table(i,5);
    end
    %====================================================================================
    
    
    % developed the stiffness matrix for an element
    for i=1:10
       k(1,1,i)= element_table(i,6)^2;
       k(1,2,i)= element_table(i,6)*element_table(i,7);
       k(1,3,i)= -1*element_table(i,6)^2;
       k(1,4,i)= -1*element_table(i,6)*element_table(i,7);
       k(2,1,i)= element_table(i,6)*element_table(i,7);
       k(2,2,i)= element_table(i,7)^2;
       k(2,3,i)= -1*element_table(i,6)*element_table(i,7);
       k(2,4,i)= -1*element_table(i,7)^2;
       k(3,1,i)= -1*element_table(i,6)^2;
       k(3,2,i)= -1*element_table(i,6)*element_table(i,7);
       k(3,3,i)= element_table(i,6)^2;
       k(3,4,i)= element_table(i,6)*element_table(i,7);
       k(4,1,i)= -1*element_table(i,6)*element_table(i,7);
       k(4,2,i)= -1*element_table(i,7)^2;
       k(4,3,i)= element_table(i,6)*element_table(i,7);
       k(4,4,i)= element_table(i,7)^2;
       k(:,:,i) = k(:,:,i).*((element_table(i,3)*element_table(i,4))/element_table(i,5));
    end

    %add the degree of freedom for each element stiffness matrix into the same degree of freedom in the structural matrix.
    K=zeros(12);
    for i=1:10
        K(2*element_table(i,1)-1,2*element_table(i,1)-1) = K(2*element_table(i,1)-1,2*element_table(i,1)-1) + k(1,1,i);
        K(2*element_table(i,1)-1,2*element_table(i,1)) = K(2*element_table(i,1)-1,2*element_table(i,1)) + k(1,2,i);
        K(2*element_table(i,1)-1,2*element_table(i,2)-1) = K(2*element_table(i,1)-1,2*element_table(i,2)-1) + k(1,3,i);
        K(2*element_table(i,1)-1,2*element_table(i,2)) = K(2*element_table(i,1)-1,2*element_table(i,2)) + k(1,4,i);
        K(2*element_table(i,1),2*element_table(i,1)-1) = K(2*element_table(i,1),2*element_table(i,1)-1) + k(2,1,i);
        K(2*element_table(i,1),2*element_table(i,1)) = K(2*element_table(i,1),2*element_table(i,1)) + k(2,2,i);
        K(2*element_table(i,1),2*element_table(i,2)-1) = K(2*element_table(i,1),2*element_table(i,2)-1) + k(2,3,i);
        K(2*element_table(i,1),2*element_table(i,2)) = K(2*element_table(i,1),2*element_table(i,2)) + k(2,4,i);   
        K(2*element_table(i,2)-1,2*element_table(i,1)-1) = K(2*element_table(i,2)-1,2*element_table(i,1)-1) + k(3,1,i);
        K(2*element_table(i,2)-1,2*element_table(i,1)) = K(2*element_table(i,2)-1,2*element_table(i,1)) + k(3,2,i);
        K(2*element_table(i,2)-1,2*element_table(i,2)-1) = K(2*element_table(i,2)-1,2*element_table(i,2)-1) + k(3,3,i);
        K(2*element_table(i,2)-1,2*element_table(i,2)) = K(2*element_table(i,2)-1,2*element_table(i,2)) + k(3,4,i);     
        K(2*element_table(i,2),2*element_table(i,1)-1) = K(2*element_table(i,2),2*element_table(i,1)-1) + k(4,1,i);
        K(2*element_table(i,2),2*element_table(i,1)) = K(2*element_table(i,2),2*element_table(i,1)) + k(4,2,i);
        K(2*element_table(i,2),2*element_table(i,2)-1) = K(2*element_table(i,2),2*element_table(i,2)-1) + k(4,3,i);
        K(2*element_table(i,2),2*element_table(i,2)) = K(2*element_table(i,2),2*element_table(i,2)) + k(4,4,i);    
    end
    K_reduced = K(1:8,1:8);  

    %the displacement of each node (Q)
    Q = K_reduced\F;
    for i=9:12
        Q(i)=0;
    end

    %the stress in each element. 
    for i=1:10
        sigma(i) = (element_table(i,3)/element_table(i,5))*[-1*element_table(i,6) -1*element_table(i,7) element_table(i,6) element_table(i,7)]*[Q(2*element_table(i,1)-1);Q(2*element_table(i,1));Q(2*element_table(i,2)-1);Q(2*element_table(i,2))];
    end
    
    K_reaction = K(9:12,1:12);  
    R = K_reaction*Q;
end