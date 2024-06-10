function [Z,E] = solve_spnlrr(Method,fun,X,lambda)
% obtain coefficients matrix
% switch Method{fun}
%
%     %--------------------Schatten-p norm minimization based dLRR by ADM--------------------------
%     case 'spdualNulrr'
%         [err_lrr1,Z,E] = lrr(X,lambda);
%          save err_mlrr1  err_lrr1;
%         % plot(err_lrr1);
%
%     case 'spdual23lrr'
%         [err_lrr2, Z,E] = nlrr_dualS23(X,lambda);
%         save err_mlrr23  err_lrr2;
%         %  plot(err_lrr2);
%
%     case 'spdual12lrr'
%         [err_lrr3,Z,E] = nlrr_dualS12(X,lambda);
%          save err_mlrr12  err_lrr3;
%         %  plot(err_lrr3);
% end

for fun = 1:3
    if fun == 1
        [err_lrr1,Z,E] = lrr(X,0.2);
        save err_mlrr1  err_lrr1;
        
    elseif fun ==2
        [err_lrr2, Z,E] = nlrr_dualS23(X,0.7);
        save err_mlrr23  err_lrr2;
        
    else
        [err_lrr3,Z,E] = nlrr_dualS12(X,0.8);
        save err_mlrr12  err_lrr3;
    end
end

plot(log(err_lrr1), 'g--', 'linewidth', 3);
hold on;
plot(log(err_lrr2), 'b-.', 'linewidth', 3);
hold on;
plot(log(err_lrr3), 'r-', 'linewidth', 3);
grid on;
legend(['p=1.0'], ['p=2/3'], ['p=1/2']);
xlabel('Number if Iterations', 'fontsize', 8);
ylabel('Error Values', 'fontsize', 8);

end




