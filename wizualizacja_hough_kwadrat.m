clc;
clear;
close all;

% Parametry kwadratu
image_size = 300; % Rozmiar obrazu
square_size = 100; % Rozmiar kwadratu
offset = 50; % Odstęp od brzegu obrazu

% Tworzenie obrazu z kwadratem
image = zeros(image_size, image_size);
image(offset:offset+square_size, offset) = 1; % Lewa krawędź
image(offset:offset+square_size, offset+square_size) = 1; % Prawa krawędź
image(offset, offset:offset+square_size) = 1; % Górna krawędź
image(offset+square_size, offset:offset+square_size) = 1; % Dolna krawędź

% Wyświetlenie oryginalnego obrazu
figure;
subplot(1, 2, 1);
imshow(image);
title('Oryginalny obraz z kwadratem');

% Transformata Hougha
[rows, cols] = size(image);
theta = linspace(-180, 179, 360); % Kąty od -180 do 179 stopni
rho_max = round(sqrt(rows^2 + cols^2));
rho = -rho_max:rho_max;
accumulator = zeros(length(rho), length(theta));

% Przeszukiwanie wszystkich pikseli
for x = 1:cols
    for y = 1:rows
        if image(y, x) == 1 % Jeśli piksel jest na krawędzi
            for t_idx = 1:length(theta)
                th = theta(t_idx) * pi / 180; % Konwersja stopni na radiany
                r = round(x * cos(th) + y * sin(th));
                rho_idx = r + rho_max + 1;
                accumulator(rho_idx, t_idx) = accumulator(rho_idx, t_idx) + 1;
            end
        end
    end
end

% Wyświetlenie akumulatora
subplot(1, 2, 2);
imshow(imadjust(rescale(accumulator)), 'XData', theta, 'YData', rho, ...
    'InitialMagnification', 'fit');
title('Akumulator Hougha');
xlabel('\theta (degrees)');
ylabel('\rho');
axis on;
axis normal;
colormap(gca, hot);

disp('Transformacja Hougha zakończona.');
