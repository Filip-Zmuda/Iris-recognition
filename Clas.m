clear all;
close all;
clc;

% Ścieżka do folderu z danymi
folder_path = 'OFTA/OFTA_reshape/codes';

% % Ustawienie ręcznego progu decyzyjnego
% threshold = 0.36;
% 
% % Inicjalizacja zmiennych
% num_folders = 21; % Liczba folderów (od 1 do 21)
% num_codes_per_eye = 3; % Liczba plików .mat na oko (o_sr11.mat, o_sr12.mat, o_sr13.mat)
% 
% all_iris_codes = cell(num_folders, num_codes_per_eye);
% all_files = cell(num_folders, num_codes_per_eye);
% 
% % Wczytanie wszystkich kodów tęczówek do struktur
% for i = 1:num_folders
%     folder_name = sprintf('%d', i);
%     files = dir(fullfile(folder_path, folder_name, '*.mat'));
%     for j = 1:num_codes_per_eye
%         file_name = files(j).name;
%         all_files{i, j} = file_name;
%         file_path = fullfile(folder_path, folder_name, file_name);
%         all_iris_codes{i, j} = load_code(file_path);
%     end
% end
% 
% correct_acceptances = 0;
% false_acceptances = 0;
% correct_rejections = 0;
% false_rejections = 0;
% authentic_comparisons = 0;
% authentic_distances = [];
% impostor_comparisons = 0;
% impostor_distances = [];
% 
% % Porównanie kodów tęczówek tego samego oka
% for i = 1:num_folders
%     for j = 1:num_codes_per_eye
%         for k = j+1:num_codes_per_eye
%             hamming_distance = calculate_hamming_distance(all_iris_codes{i, j}, all_iris_codes{i, k});
%             classification = classify_code(hamming_distance, threshold);
% 
%             % Sprawdzenie poprawności klasyfikacji na podstawie nazw plików
%             file1 = all_files{i, j};
%             file2 = all_files{i, k};
% 
%             if is_same_eye(file1, file2)
%                 authentic_comparisons = authentic_comparisons + 1;
%                 authentic_distances = [authentic_distances, hamming_distance];
%                 if strcmp(classification, 'AKCEPTACJA')
%                     correct_acceptances = correct_acceptances + 1;
%                 else
%                     false_rejections = false_rejections + 1;
%                 end
%             else
%                 if strcmp(classification, 'AKCEPTACJA')
%                     false_acceptances = false_acceptances + 1;
%                 else
%                     correct_rejections = correct_rejections + 1;
%                 end
%             end
%         end
%     end
% end
% 
% if authentic_comparisons > 0
%     FRR = false_rejections / authentic_comparisons;
% else
%     FRR = 0;
% end
% 
% % Porównanie kodów tęczówek różnych oczu
% for i = 1:num_folders
%     for j = 1:num_codes_per_eye
%         for k = i+1:num_folders
%             for l = 1:num_codes_per_eye
%                 hamming_distance = calculate_hamming_distance(all_iris_codes{i, j}, all_iris_codes{k, l});
%                 classification = classify_code(hamming_distance, threshold);
% 
%                 % Sprawdzenie poprawności klasyfikacji na podstawie nazw plików
%                 file1 = all_files{i, j};
%                 file2 = all_files{k, l};
% 
%                 if is_same_eye(file1, file2)
%                     if strcmp(classification, 'AKCEPTACJA')
%                         correct_acceptances = correct_acceptances + 1;
%                     else
%                         false_rejections = false_rejections + 1;
%                     end
%                 else
%                     impostor_comparisons = impostor_comparisons + 1;
%                     impostor_distances = [impostor_distances, hamming_distance];
%                     if strcmp(classification, 'AKCEPTACJA')
%                         false_acceptances = false_acceptances + 1;
%                     else
%                         correct_rejections = correct_rejections + 1;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% if impostor_comparisons > 0
%     FAR = false_acceptances / impostor_comparisons;
% else
%     FAR = 0;
% end
% 
% % Wyświetlanie wyników
% fprintf('Liczba folderów: %d\n', num_folders);
% fprintf('Próg decyzyjny: %.2f\n', threshold);
% fprintf('Poprawne akceptacje (CA): %d\n', correct_acceptances);
% fprintf('Fałszywe akceptacje (FA): %d\n', false_acceptances);
% fprintf('Poprawne odrzucenia (CR): %d\n', correct_rejections);
% fprintf('Fałszywe odrzucenia (FR): %d\n', false_rejections);
% fprintf('Autentyczne porównania: %d\n', authentic_comparisons);
% fprintf('Porównania impostorów: %d\n', impostor_comparisons);
% 
% mean_authentic = mean(authentic_distances);
% var_authentic = var(authentic_distances);
% mean_impostor = mean(impostor_distances);
% var_impostor = var(impostor_distances);
% 
% % Wyświetlenie średnich i wariancji
% disp(['Średnia (autentyczne): ', num2str(mean_authentic)]);
% disp(['Wariancja (autentyczne): ', num2str(var_authentic)]);
% disp(['Średnia (impostorzy): ', num2str(mean_impostor)]);
% disp(['Wariancja (impostorzy): ', num2str(var_impostor)]);
% 
% % Generowanie funkcji gęstości Gaussa
% x = linspace(0, 1, 1000);
% pdf_authentic = normpdf(x, mean_authentic, sqrt(var_authentic));
% pdf_impostor = normpdf(x, mean_impostor, sqrt(var_impostor));
% 
% % Rysowanie funkcji gęstości Gaussa
% figure;
% hold on;
% plot(x, pdf_authentic, 'b', 'LineWidth', 2);
% plot(x, pdf_impostor, 'r', 'LineWidth', 2);
% fill([x fliplr(x)], [pdf_authentic zeros(size(pdf_authentic))], 'b', 'FaceAlpha', 0.3);
% fill([x fliplr(x)], [pdf_impostor zeros(size(pdf_impostor))], 'r', 'FaceAlpha', 0.3);
% xlabel('Odległość');
% ylabel('Prawdopodobieństwo');
% legend('P_{u}(x) (autentyczne)', 'P_{o}(x) (impostorzy)');
% title('Funkcje gęstości prawdopodobieństwa Gaussa');
% grid on;
% hold off;
% 
% % Wyświetlenie FAR i FRR
% fprintf('FAR: %.6f\n', FAR);
% fprintf('FRR: %.6f\n', FRR);
% 
% % Wykres DET
% figure;
% plot(FAR, FRR, '-o');
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% xlabel('FAR');
% ylabel('FRR');
% title('DET - Detection Error Trade-off');
% grid on;







% Inicjalizacja zmiennych
num_folders = 21; % Liczba folderów (od 1 do 21)
num_codes_per_eye = 3; % Liczba plików .mat na oko (o_sr11.mat, o_sr12.mat, o_sr13.mat)

% Utworzenie wektora z różnymi wartościami progu decyzyjnego
thresholds = 0:0.01:0.5;

% Inicjalizacja wektorów FAR i FRR
FAR = zeros(size(thresholds));
FRR = zeros(size(thresholds));

% Wczytanie wszystkich kodów tęczówek do struktur
all_iris_codes = cell(num_folders, num_codes_per_eye);
all_files = cell(num_folders, num_codes_per_eye);
for i = 1:num_folders
    folder_name = sprintf('%d', i);
    files = dir(fullfile(folder_path, folder_name, '*.mat'));
    for j = 1:num_codes_per_eye
        file_name = files(j).name;
        all_files{i, j} = file_name;
        file_path = fullfile(folder_path, folder_name, file_name);
        all_iris_codes{i, j} = load_code(file_path);
    end
end

% Obliczenie FAR i FRR dla każdego progu decyzyjnego
for t = 1:length(thresholds)
    threshold = thresholds(t);
    
    correct_acceptances = 0;
    false_acceptances = 0;
    correct_rejections = 0;
    false_rejections = 0;
    authentic_comparisons = 0;
    authentic_distances = [];
    impostor_comparisons = 0;
    impostor_distances = [];
    
    % Porównanie kodów tęczówek tego samego oka
    for i = 1:num_folders
        for j = 1:num_codes_per_eye
            for k = j+1:num_codes_per_eye
                hamming_distance = calculate_hamming_distance(all_iris_codes{i, j}, all_iris_codes{i, k});
                classification = classify_code(hamming_distance, threshold);
                
                % Sprawdzenie poprawności klasyfikacji na podstawie nazw plików
                file1 = all_files{i, j};
                file2 = all_files{i, k};
                
                if is_same_eye(file1, file2)
                    authentic_comparisons = authentic_comparisons + 1;
                    authentic_distances = [authentic_distances, hamming_distance];
                    if strcmp(classification, 'AKCEPTACJA')
                        correct_acceptances = correct_acceptances + 1;
                    else
                        false_rejections = false_rejections + 1;
                    end
                else
                    if strcmp(classification, 'AKCEPTACJA')
                        false_acceptances = false_acceptances + 1;
                    else
                        correct_rejections = correct_rejections + 1;
                    end
                end
            end
        end
    end
    
    if authentic_comparisons > 0
        FRR(t) = false_rejections / authentic_comparisons;
    else
        FRR(t) = 0;
    end
    
    % Porównanie kodów tęczówek różnych oczu
    for i = 1:num_folders
        for j = 1:num_codes_per_eye
            for k = i+1:num_folders
                for l = 1:num_codes_per_eye
                    hamming_distance = calculate_hamming_distance(all_iris_codes{i, j}, all_iris_codes{k, l});
                    classification = classify_code(hamming_distance, threshold);
                    
                    % Sprawdzenie poprawności klasyfikacji na podstawie nazw plików
                    file1 = all_files{i, j};
                    file2 = all_files{k, l};
                    
                    if is_same_eye(file1, file2)
                        if strcmp(classification, 'AKCEPTACJA')
                            correct_acceptances = correct_acceptances + 1;
                        else
                            false_rejections = false_rejections + 1;
                        end
                    else
                        impostor_comparisons = impostor_comparisons + 1;
                        impostor_distances = [impostor_distances, hamming_distance];
                        if strcmp(classification, 'AKCEPTACJA')
                            false_acceptances = false_acceptances + 1;
                        else
                            correct_rejections = correct_rejections + 1;
                        end
                    end
                end
            end
        end
    end
    
    if impostor_comparisons > 0
        FAR(t) = false_acceptances / impostor_comparisons;
    else
        FAR(t) = 0;
    end
end

% Rysowanie krzywej DET
figure;
plot(FAR, FRR, '-o');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('FAR');
ylabel('FRR');
title('DET - Detection Error Trade-off');
grid on;



function code = load_code(file_path)
    % Wczytuje kod tęczówki z pliku .mat
    data = load(file_path);
    % Zakładamy, że zmienna w pliku .mat zawiera ciąg zer i jedynek jako logical(1x160896)
    % Sprawdź nazwę zmiennej w pliku .mat
    var_names = fieldnames(data);
    if numel(var_names) == 1
        % Wczytaj pierwszą zmienną w pliku .mat
        code = data.(var_names{1});
    else
        error('Nieprawidłowa liczba zmiennych w pliku .mat.');
    end
end


function hamming_distance = calculate_hamming_distance(code1, code2)
    % Oblicza odległość Hamminga między dwoma kodami tęczówek
    hamming_distance = sum(code1 ~= code2) / length(code1);
end

function classification = classify_code(hamming_distance, threshold)
    % Klasyfikuje kod tęczówki na podstawie odległości Hamminga i zadanego progu
    if hamming_distance <= threshold
        classification = 'AKCEPTACJA';
    else
        classification = 'ODRZUCENIE';
    end
end

function result = is_same_eye(file1, file2)
    % Sprawdza, czy dwa pliki pochodzą od tego samego oka
    % Zakładamy, że nazwa pliku ma format SXXXXXLXX
    prefix1 = file1(1:6);
    prefix2 = file2(1:6);
    result = strcmp(prefix1, prefix2);
end
