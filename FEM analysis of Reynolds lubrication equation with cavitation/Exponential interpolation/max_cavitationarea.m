%determine max abs of cavitation points. These errors are saved in
%max_contatc.mat.

% n = [3:6];
% for j = 1:length(n)
%     nRef = n(j);
%     LABO;
%     num_iter(j) = ii; % num_iter = [10    12    15    21]
%     p3 = full(p3);
%     neg = nonzeros(setdiff(p3,p3(xx)));
%     e(j) = max(abs(neg));
%     
% end
%%
% plot errors fitted with an exponential interpolation
e = load ('max_contact.mat');
e = e.e;
n = [3:6];
logp1 = polyfit([n(1):n(end)],log10(e),1);
logpred1 = 10.^polyval(logp1,[n(1):0.01:n(end)]);

plot([n(1):0.01:n(end)],logpred1,'LineWidth',1)
hold on
plot([n(1):n(end)],e,'*','LineWidth',1)
grid on
xlabel('nRef')
ylabel('max|error|')
text(n(1),e(1),num2str(e(1)))
text(n(2),e(2),num2str(e(2)))
text(n(3),e(3),num2str(e(3)))
text(n(4),e(4),num2str(e(4)))