function f = plot_fun(r1,r2,L,density)
    node_coor = [2*L L; 2*L 0; L L; L 0; 0 L; 0 0];
    element_table = [3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4];
    f=0;
    for i=1:10
        if i<=6
            area = pi*r1*r1;
        end

        if i>6
            area = pi*r2*r2;
        end

        x_diff = (node_coor(element_table(i,2),1)- node_coor(element_table(i,1),1));
        y_diff = (node_coor(element_table(i,2),2)- node_coor(element_table(i,1),2));
        length = sqrt(x_diff^2 + y_diff^2);
        f = f + area*length*density;
    end
end