function iris_code = iris_recognition_ofta(image_path)
%% Wczytywanie, preprocessing, binaryzacja

oko_org=imread(image_path);

oko_org=rgb2gray(oko_org);
dzielnik=2;

%oko_resized=imresize(oko_org,[480/dzielnik,640/dzielnik], 'bilinear');     %csil
oko_org=imresize(oko_org,[600,860], 'bilinear');
oko_resized=imresize(oko_org,[600/dzielnik,860/dzielnik], 'bilinear');      %ofta

oko_med = median_filter(oko_resized, 3);
oko_med = median_filter(oko_med, 3);
oko_blurred=custom_gaussian_filter(oko_med, 1.8);

%% binaryzacja

oko=double(oko_blurred);

sx= [-1 0 1;
     -2 0 2;
     -1 0 1];
sy= [-1 -2 -1;
     0 0 0;
     1 2 1];

gradient_x = zeros(size(oko) - [2, 2]);
gradient_y = zeros(size(oko) - [2, 2]);

for i = 1:size(oko, 1) - 2 
    for j = 1:size(oko, 2) - 2 
        Gx = sum(sum(sx.*oko(i:i+2, j:j+2))); 
        Gy = sum(sum(sy.*oko(i:i+2, j:j+2))); 
        gradient_x(i+1, j+1) = Gx;
        gradient_y(i+1, j+1) = Gy;
        grad(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
    end
end

%treshold = 64/dzielnik;
treshold = 32/dzielnik;
grad = uint8(grad);
grad=grad>treshold;

%% Transformata Hougha
[rows, cols] = size(grad);

% Inicjalizacja parametrów
% min_iris_radius = round(90/dzielnik); % Minimalny promień tęczówki
% max_iris_radius = round(140/dzielnik); % Maksymalny promień tęczówki
% min_pupil_radius = round(50/dzielnik); % Minimalny promień źrenicy
% max_pupil_radius = round(70/dzielnik); % Maksymalny promień źrenicy     %csil

min_pupil_radius = round(60/dzielnik); % Minimalny promień źrenicy
max_pupil_radius = round(110/dzielnik); % Maksymalny promień źrenicy
min_iris_radius = round((height(oko_resized)/2)-20); % Minimalny promień tęczówki
max_iris_radius = round((height(oko_resized)/2)+10); % Maksymalny promień tęczówki     %%ofta

% Inicjalizacja macierzy akumulatorów Hougha dla tęczówki i źrenicy
iris_accumulator = zeros(rows, cols, max_iris_radius);
pupil_accumulator = zeros(rows, cols, max_pupil_radius);

tic;

% Detekcja tęczówki
for r = min_iris_radius:max_iris_radius % Dla każdego promienia w zakresie dla tęczówki
    for x = 1:cols
        for y = 1:rows
            if grad(y, x) == 1 % Jeśli piksel jest na krawędzi (gradient jest równy 1)
                for theta = 0:pi/50:2*pi % Dla każdego kąta theta
                    a = round(x - r * cos(theta)); % Oblicz współrzędną x środka okręgu
                    b = round(y - r * sin(theta)); % Oblicz współrzędną y środka okręgu
                    if a > 0 && a <= cols && b > 0 && b <= rows % Sprawdź, czy współrzędne są w granicach obrazu
                        iris_accumulator(b, a, r) = iris_accumulator(b, a, r) + 1; % Inkrementuj wartość akumulatora Hougha dla tęczówki
                    end
                end
            end
        end
    end
end

% Znajdź najjaśniejszy punkt w akumulatorze dla tęczówki
[max_iris_votes, max_iris_index] = max(iris_accumulator(:));
[iris_y, iris_x, iris_r] = ind2sub(size(iris_accumulator), max_iris_index);

% Detekcja źrenicy (w ograniczonym obszarze)
for r = min_pupil_radius:max_pupil_radius % Dla każdego promienia w zakresie dla źrenicy
    for x = 1:cols
        for y = 1:rows
            if grad(y, x) == 1 % Jeśli piksel jest na krawędzi (gradient jest równy 1)
                for theta = 0:pi/50:2*pi % Dla każdego kąta theta
                    a = round(x - r * cos(theta)); % Oblicz współrzędną x środka okręgu
                    b = round(y - r * sin(theta)); % Oblicz współrzędną y środka okręgu
                    if a > 0 && a <= cols && b > 0 && b <= rows % Sprawdź, czy współrzędne są w granicach obrazu
                        pupil_accumulator(b, a, r) = pupil_accumulator(b, a, r) + 1; % Inkrementuj wartość akumulatora Hougha dla źrenicy
                    end
                end
            end
        end
    end
end

% Znajdź najjaśniejszy punkt w akumulatorze dla źrenicy
[max_pupil_votes, max_pupil_index] = max(pupil_accumulator(:));
[pupil_y, pupil_x, pupil_r] = ind2sub(size(pupil_accumulator), max_pupil_index);

elapsed_time=toc;
fprintf('Czas wykonania: %.4f sekundy\n', elapsed_time);


%% Algorytm Daugmana

% % Parametry transformacji Daugmana
% [rows, cols] = size(grad);
% % rmin = round(50/dzielnik); % Minimalny promień źrenicy
% % rmax = round(70/dzielnik); % Maksymalny promień źrenicy
% % irisrmin = round(90/dzielnik); % Minimalny promień tęczówki
% % irisrmax = round(140/dzielnik); % Maksymalny promień tęczówki %%csil
% 
% rmin = round(50/dzielnik); % Minimalny promień źrenicy
% rmax = round(110/dzielnik); % Maksymalny promień źrenicy
% irisrmin = round((height(oko_resized)/2)-80); % Minimalny promień tęczówki
% irisrmax = round((height(oko_resized)/2)+10); % Maksymalny promień tęczówki     %%ofta
% 
% pupil_accumulator = zeros(rows, cols, rmax - rmin + 1);
% iris_accumulator = zeros(rows, cols, irisrmax - irisrmin + 1);
% 
% tic;
% 
% %Detekcja tęczówki
% for x = 1:cols
%     for y = 1:rows
%         for r = irisrmin:irisrmax
%             for theta = 0:pi/6:2*pi
%                 a = round(x - r * cos(theta));
%                 b = round(y - r * sin(theta));
%                 if a > 0 && a <= cols && b > 0 && b <= rows
%                     iris_accumulator(y, x, r - irisrmin + 1) = iris_accumulator(y, x, r - irisrmin + 1) + grad(b, a);
%                 end
%             end
%         end
%     end
% end
% 
% % Znajdź najjaśniejszy punkt w akumulatorze dla tęczówki
% [max_iris_votes, max_iris_index] = max(iris_accumulator(:));
% [iris_y, iris_x, iris_r] = ind2sub(size(iris_accumulator), max_iris_index);
% iris_r = iris_r + irisrmin - 1;
% 
% % Detekcja źrenicy
% % for r = rmin:rmax
% % for x = (iris_x - iris_r + r):(iris_x + iris_r - r)
% %     for y = (iris_y - iris_r + r):(iris_y + iris_r - r)
% %        if (x - iris_x)^2 + (y - iris_y)^2 <= (iris_r^2) % Sprawdź, czy punkt (x, y) jest wewnątrz okręgu tęczówki
% % 
% %             for theta = 0:pi/18:2*pi
% %                 a = round(x - r * cos(theta));
% %                 b = round(y - r * sin(theta));
% %                 if a > 0 && a <= cols && b > 0 && b <= rows
% %                     pupil_accumulator(y, x, r - rmin + 1) = pupil_accumulator(y, x, r - rmin + 1) + grad(b, a);
% %                 end
% %             end
% %         end
% %     end
% %     end
% % end
% 
% % Detekcja źrenicy
% for x = 1:cols
%     for y = 1:rows
%         for r = rmin:rmax
%             for theta = 0:pi/6:2*pi
%                 a = round(x - r * cos(theta));
%                 b = round(y - r * sin(theta));
%                 if a > 0 && a <= cols && b > 0 && b <= rows
%                     pupil_accumulator(y, x, r - rmin + 1) = pupil_accumulator(y, x, r - rmin + 1) + grad(b, a);
%                 end
%             end
%         end
%     end
% end
% 
% % Znajdź najjaśniejszy punkt w akumulatorze dla źrenicy
% [max_pupil_votes, max_pupil_index] = max(pupil_accumulator(:));
% [pupil_y, pupil_x, pupil_r] = ind2sub(size(pupil_accumulator), max_pupil_index);
% pupil_r = pupil_r + rmin - 1;
% 
% 
% 
% elapsed_time=toc;
% fprintf('Czas wykonania: %.4f sekundy\n', elapsed_time);

%% Rysowanie okręgów

figure(1);
imshow(grad);
hold on;
iris_theta = 0:0.01:(2*pi);
pupil_theta = 0:0.01:(2*pi);
plot(iris_x + iris_r * cos(iris_theta), iris_y + iris_r * sin(iris_theta), 'r', 'LineWidth', 2);
plot(pupil_x + pupil_r * cos(pupil_theta), pupil_y + pupil_r * sin(pupil_theta), 'g', 'LineWidth', 2);
title('Znalezione tęczówka i źrenica');
hold off;

figure(2);
imshow(oko_org);
hold on;
plot(dzielnik*iris_x + dzielnik*iris_r * cos(iris_theta), dzielnik*iris_y + dzielnik*iris_r * sin(iris_theta), 'r', 'LineWidth', 2);
plot(dzielnik*pupil_x + dzielnik*pupil_r * cos(pupil_theta), dzielnik*pupil_y + dzielnik*pupil_r * sin(pupil_theta), 'g', 'LineWidth', 2);
title('Znalezione tęczówka i źrenica');
hold off;

%% Normalizacja tęczówki
theta = 0:0.005:2*pi;
radii = linspace(pupil_r*dzielnik, iris_r*dzielnik, 128);
[thetaGrid, radiiGrid] = meshgrid(theta, radii);
% Konwersja współrzędnych biegunowych do kartezjańskich
x = pupil_x*dzielnik + radiiGrid .* cos(thetaGrid);
y = pupil_y*dzielnik + radiiGrid .* sin(thetaGrid);

normalizedIris = interp2(double(oko_org), x, y, 'linear', 0);

%% Filtr Gabora

% Parametry filtra Gabora
wavelength = 8;      % Długość fali
theta = 0;    % Orientacja w radianach
sigmaX = 4;          % Standardowe odchylenie w osi X
sigmaY = 4;          % Standardowe odchylenie w osi Y
phaseOffset = 0;     % Przesunięcie fazowe

% Definicja siatki współrzędnych
[x, y] = meshgrid(-fix(sigmaX*3):fix(sigmaX*3), -fix(sigmaY*3):fix(sigmaY*3)); %generacja


% Obrót współrzędnych
xPrime = x .* cos(theta) + y .* sin(theta);
yPrime = -x .* sin(theta) + y .* cos(theta);

% Generowanie rzeczywistej i urojonej części filtra
realPart = exp(-0.5 * (xPrime.^2 / sigmaX^2 + yPrime.^2 / sigmaY^2)) .* ...
               cos(2 * pi * xPrime / wavelength + phaseOffset);
imagPart = exp(-0.5 * (xPrime.^2 / sigmaX^2 + yPrime.^2 / sigmaY^2)) .* ...
               sin(2 * pi * xPrime / wavelength + phaseOffset);

% Tworzenie filtra zespolonego
gaborFilter = realPart + 1i * imagPart;

realFilteredImg = filter2(real(gaborFilter), normalizedIris);
imagFilteredImg = filter2(imag(gaborFilter), normalizedIris);
filteredImgComplex = realFilteredImg + 1i * imagFilteredImg;

    %% IrisCode

    % Binaryzacja
    threshold = mean(filteredImgComplex(:));
    binary_image = filteredImgComplex > threshold;

    % Kodowanie zbinaryzowanego obrazu do jednowymiarowego wektora
    iris_code = reshape(binary_image, 1, []);

    % Konwersja do double a następnie do string
    %iris_code_str = sprintf('%d', double(iris_code));
    iris_code_image = reshape((iris_code.'), size(filteredImgComplex));

    % Sprawdzanie długości kodu
    expected_length = numel(filteredImgComplex);
    actual_length = length(iris_code);
    if expected_length == actual_length
        disp('The length of the iris code is correct.');
    else
        disp('The length of the iris code is incorrect.');
    end

%     binary_image_imag = realFilteredImg > 0;
% binary_image_real = imagFilteredImg > 0;
% binary_image = [binary_image_real;binary_image_imag];
% % Encode the binary image into a bit string
% iris_code = reshape(binary_image, 1, []);
% 
% iris_code_image = reshape( iris_code,128,[]);

figure(3);
imshow(iris_code_image, []);
title('Iris Code Image');
end
