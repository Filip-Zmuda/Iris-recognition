clear all;
close all;
clc;

% Ścieżka do folderu z kodami tęczówek
codes_path = 'OFTA\OFTA_reshape\codes';
threshold = 0.38; % Próg decyzyjny dla odległości Hamminga

% Inicjalizacja macierzy pomyłek
num_people = 21;
confusion_matrix = zeros(num_people, num_people);

% Inicjalizacja liczników dla FAR i FRR
total_comparisons = 0;
false_accepts = 0;
false_rejects = 0;

% Inicjalizacja zmiennych do przechowywania wyników odległości Hamminga
hamming_distances = cell(num_people, 3); % 3 zdjęcia dla każdej osoby

tic;

% Przetwarzanie bazy danych
for i = 1:num_people
    person_folder = fullfile(codes_path, num2str(i));
    code_files = dir(fullfile(person_folder, '*.mat'));
    
    for j = 1:length(code_files)
        % Załaduj kod tęczówki z bazy danych
        load(fullfile(person_folder, code_files(j).name), 'iris_code');
        
        % Zapisz kod tęczówki do macierzy odległości Hamminga
        hamming_distances{i, j} = iris_code;
    end
end

% Obliczenie odległości Hamminga i ocena systemu
for i = 1:num_people
    for j = 1:3
        new_iris_code = hamming_distances{i, j};
        
        % Obliczenie odległości Hamminga dla każdego kodu tęczówki w bazie danych
        for k = 1:num_people
            if k ~= i
                for l = 1:3
                    reference_iris_code = hamming_distances{k, l};

                    if ~isempty(new_iris_code) && ~isempty(reference_iris_code)
                        % Oblicz odległość Hamminga
                        hammingDist = hammingDistance(new_iris_code, reference_iris_code);
                        total_comparisons = total_comparisons + 1;

                        % Aktualizacja macierzy pomyłek
                        if hammingDist < threshold
                            confusion_matrix(i, k) = confusion_matrix(i, k) + 1;
                            false_accepts = false_accepts + 1;
                        end
                    end
                end
            end
        end

        % Sprawdzenie samego siebie
        for l = 1:3
            if l ~= j
                reference_iris_code = hamming_distances{i, l};

                if ~isempty(new_iris_code) && ~isempty(reference_iris_code)
                    % Oblicz odległość Hamminga
                    hammingDist = hammingDistance(new_iris_code, reference_iris_code);

                    % Aktualizacja macierzy pomyłek
                    if hammingDist < threshold
                        confusion_matrix(i, i) = confusion_matrix(i, i) + 1;
                    else
                        false_rejects = false_rejects + 1;
                    end
                end
            end
        end
    end
end

elapsed_time = toc;
fprintf('Czas szukania: %.4f sekundy\n', elapsed_time);

% Obliczenie FAR, FRR
total_accepts = sum(confusion_matrix(:));
FAR = false_accepts / total_comparisons;
FRR = false_rejects / total_accepts;

% Wyświetlenie wyników
fprintf('FAR: %.4f\n', FAR);
fprintf('FRR: %.4f\n', FRR);

% Wyświetlenie macierzy pomyłek
figure;
imagesc(confusion_matrix); % Wizualizacja macierzy pomyłek
colorbar; % Dodanie paska kolorów do wskazania zakresu wartości
xlabel('Przewidywana osoba'); % Opis osi X
ylabel('Rzeczywista osoba'); % Opis osi Y
title('Macierz pomyłek');
