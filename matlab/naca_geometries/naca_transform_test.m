function [] = naca_transform_test()

Set_Default_Plot_Properties();

m = 12/100; % Max camber (x/100)
p = 5/10;   % Location of max camber (y/10)
t = 15/100; % Max thickness over chord (zz / 100)
c = 1.0;    % Chord length
a = 8;     % Angle of attack

x0 = [linspace(0.00001,0.1,20), linspace(0.11,1,80)];
[~, xu0, yu0, xl0, yl0, yc] = naca_coords(x0, 0.0, 0.0, 0.12, 0.0);
% [~, xu0, yu0, xl0, yl0, yc] = naca_coords(x0, 0.04, 0.4, 0.12, 0.0);

xu = nan(size(x0));
yu = nan(size(x0));
xl = nan(size(x0));
yl = nan(size(x0));
for i = 1:length(x0)
    % Upper surface.
    x = xu0(i);
    y = yu0(i);
    [dx, dy] = naca_transform(x, y, m, p, t, c, a);
    xu(i) = x + dx;
    yu(i) = y + dy;
    % Lower surface.
    x = xl0(i);
    y = yl0(i);
    [dx, dy] = naca_transform(x, y, m, p, t, c, a);
    xl(i) = x + dx;
    yl(i) = y + dy;
end

figure();
hold on;
plot([x0;xu],[yu0;yu],'k-', 'LineWidth', 1); % Upper deltas
plot([x0;xl],[yl0;yl],'k-', 'LineWidth', 1); % Lower deltas
plot(x0,yu0,'r-', 'LineWidth', 3); % Original upper
plot(x0,yl0,'r-', 'LineWidth', 3); % Original lower
% plot(x0,yc,'k-', 'LineWidth', 4); % Camber
plot(xu,yu,'b-', 'LineWidth', 3); % Transformed upper
plot(xl,yl,'b-', 'LineWidth', 3); % Transformed lower
% title(sprintf('m = %.2f, p = %.2f, t = %.2f, c = %.2f, a = %.1f',m,p,t,c,a));
xlabel('x');
ylabel('y');
axis equal;
xlim([-0.05,1.05]);
ylim([-0.2, 0.2]);

end





