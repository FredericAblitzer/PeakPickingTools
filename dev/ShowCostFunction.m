function ShowCostFunction(handles)

global s f

n = get(handles.listboxModes,'Value');

X0 = [s.OMEGA_n_0(n) 0.05];

[X1,X2] = meshgrid(linspace(0.8,1.2,51)*X0(1),logspace(-3,0,51));
RES = NaN*zeros(size(X1));


for i = 1:size(X1,1)
    i
    for j = 1:size(X1,2)
        RES(i,j) = f.FindAmplitudes(...
            [ s.OMEGA_n([1:n-1]) ; X1(i,j) ; s.OMEGA_n([n+1:end])],...
            [ s.XI_n([1:n-1]) ; X2(i,j) ; s.XI_n([n+1:end])]);
    end
end

Xsol = [s.OMEGA_n(n) s.XI_n(n)];

figure
clf
hold on
pcolor(X1/(2*pi),X2,RES)
plot(X0(1)/(2*pi),X0(2),'wo')
plot(Xsol(1)/(2*pi),Xsol(2),'mx','MarkerSize',14,'Linewidth',2)
shading interp
colormap(gray)
set(gca,'YScale','log')
xlim(X1([1 end])/(2*pi))
ylim(X2([1 end]))
xlabel('\omega_n/2\pi (Hz)')
ylabel('\xi_n')
box on

set(gca,'Layer','top')


set(gcf,'Units','pixels')
set(gcf,'Position',[0 0 500 300])

set(findall(get(gcf,'Children'),'-property','FontSize'),'Fontsize',14)

