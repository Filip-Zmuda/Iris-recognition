clear all;
close all;
clc;

% % Przykład użycia skryptu
% image_path = 'OFTA\OFTA_reshape\pictures\11\o_sr111.bmp';
% codes_path = 'OFTA\OFTA_bez_1\codes';
% threshold = 0.37; % Próg decyzyjny dla odległości Hamminga
% 
% tic;
% 
% [matched_person, minHammingDistance, access_granted] = identify_iris(image_path, codes_path, threshold);
% 
% elapsed_time=toc;
% fprintf('Czas szukania: %.4f sekundy\n', elapsed_time);

evaluate_iris_model('OFTA/OFTA_reshape/codes', 0.40)