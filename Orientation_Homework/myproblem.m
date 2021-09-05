function [X fval history] = myproblem(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,history,i)
    options = optimset('PlotFcns', @myoutput);
    [X fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
    function stop = myoutput(X,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history;i X nonlcon(X) fun(X)];
          i=i+1;
        end
    end   
end
