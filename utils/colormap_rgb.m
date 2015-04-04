function [cols] = colormap_rgb()
t = [0 139 159; 255 146 0; 47 165 80 ]./255;
R = [linspace(t(1,1),t(2,1),9) linspace(t(2,1),t(3,1),9)];
G = [linspace(t(1,2),t(2,2),9) linspace(t(2,2),t(3,2),9)];
B = [linspace(t(1,3),t(2,3),9) linspace(t(2,3),t(3,3),9)];
cols = [R([1:9 11:18])', G([1:9 11:18])', B([1:9 11:18])'];