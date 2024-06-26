clc;
clear;
close all;

% Parametry okręgu
radius = 300; % Znany promień okręgu
center_x = 900; % Współrzędna x środka okręgu
center_y = 900; % Współrzędna y środka okręgu
image_size = 1800; % Rozmiar obrazu

% Tworzenie obrazu z jednym okręgiem
image = zeros(image_size, image_size);
theta = linspace(0, 2*pi, 36000); % Kąty dla rysowania okręgu
x_circle = center_x + radius * cos(theta);
y_circle = center_y + radius * sin(theta);
for i = 1:length(x_circle)
    x = round(x_circle(i));
    y = round(y_circle(i));
    if x > 0 && x <= image_size && y > 0 && y <= image_size
        image(y, x) = 1;
    end
end


% Wyświetlenie oryginalnego obrazu
figure;
subplot(1, 2, 1);
imshow(image);
title('Oryginalny obraz z okręgiem');

% Kołowa transformata Hougha
[rows, cols] = size(image);
accumulator = zeros(rows, cols);

% Przeszukiwanie wszystkich pikseli
for x = 1:cols
    for y = 1:rows
        if image(y, x) == 1 % Jeśli piksel jest na krawędzi (należy do okręgu)
            for theta = linspace(0, 2*pi, 35) % Dla każdego kąta theta
                a = round(x - radius * cos(theta)); % Oblicz współrzędną x środka okręgu
                b = round(y - radius * sin(theta)); % Oblicz współrzędną y środka okręgu
                if a > 0 && a <= cols && b > 0 && b <= rows % Sprawdź, czy współrzędne są w granicach obrazu
                    accumulator(b, a) = accumulator(b, a) + 1; % Inkrementuj wartość akumulatora Hougha
                end
            end
        end
    end
end

% Wyświetlenie akumulatora
subplot(1, 2, 2);
imshow(imadjust(mat2gray(accumulator)));
hold on;
title('Akumulator Hougha');

% Naniesienie okręgów na akumulator
[accum_max, ind] = max(accumulator(:));
[y_peak, x_peak] = ind2sub(size(accumulator), ind);

% Rysowanie wykrytego okręgu
viscircles([x_peak, y_peak], radius, 'EdgeColor', 'b');
hold off;

disp(['Środek wykrytego okręgu: (', num2str(x_peak), ', ', num2str(y_peak), ')']);
disp(['Promień wykrytego okręgu: ', num2str(radius)]);
