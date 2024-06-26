clear all;
close all;
clc;

% iris_code1 = iris_recognition_csil('CSIL/050/R/S2050R20.jpg');
% %disp(iris_code1);
% iris_code2 = iris_recognition_csil('CSIL/050/R/S2050R10.jpg');
% %iris_code2 = iris_recognition('CSIL/030/R/S2030R18.jpg');

% iris_code1 = iris_recognition_ofta('OFTA/o_sr101.bmp');
% iris_code2 = iris_recognition_ofta('OFTA/o_sr102.bmp');
tic;
iris_code1 = iris_recognition_ofta('OFTA/o_sr81.bmp');
elapsed_time=toc;
fprintf('Czas jednego oka: %.4f sekundy\n', elapsed_time);
%iris_code2 = iris_recognition_ofta('OFTA/o_sr103.bmp');
% Ładowanie pliku .mat do struktury
dataStruct = load('OFTA/OFTA_bez_1/codes/1/o_sr12.mat');

% Pobieranie nazw pól w strukturze
fields = fieldnames(dataStruct);

% Przypisanie zawartości jedynego pola do zmiennej
iris_code2 = dataStruct.(fields{1});

%disp(iris_code2)
tic;
distance=hammingDistance(iris_code1,iris_code2);

elapsed_time=toc;
fprintf('Czas porównywania: %.4f sekundy\n', elapsed_time);

disp(distance)