function evaluate_iris_model(codes_folder, threshold)
    % codes_folder: folder z kodami tęczówek (np. 'OFTA/OFTA_reshape/codes/')
    % threshold: próg decyzyjny dla przyznania dostępu

    num_persons = 21; % Liczba osób w bazie danych
    images_per_person = 3; % Liczba zdjęć na osobę

    % Tablice do przechowywania wyników
    genuine_scores = [];
    impostor_scores = [];

    % Ocena dla każdej osoby
    for i = 1:num_persons
        for j = 1:images_per_person
            code_path = fullfile(codes_folder, num2str(i), sprintf('o_sr%d%d.mat', i, j));
            load(code_path, 'iris_code');
            new_iris_code = iris_code;
            
            for k = 1:num_persons
                for l = 1:images_per_person
                    code_path_compare = fullfile(codes_folder, num2str(k), sprintf('o_sr%d%d.mat', k, l));
                    load(code_path_compare, 'iris_code');
                    
                    if i == k && j == l
                        continue; % Pominięcie porównania tego samego kodu
                    end
                    
                    hammingDist = hammingDistance(new_iris_code, iris_code);
                    
                    if i == k
                        genuine_scores = [genuine_scores, hammingDist];
                    else
                        impostor_scores = [impostor_scores, hammingDist];
                    end
                end
            end
        end
    end

    % Obliczenie parametrów FAR i FRR
    FAR = sum(impostor_scores < threshold) / length(impostor_scores);
    FRR = sum(genuine_scores >= threshold) / length(genuine_scores);

    % Obliczenie współczynnika dyskryminacji D
    mean_genuine = mean(genuine_scores);
    mean_impostor = mean(impostor_scores);
    std_genuine = std(genuine_scores);
    std_impostor = std(impostor_scores);
    D = (mean_impostor - mean_genuine) / sqrt((std_genuine^2 + std_impostor^2) / 2);

    % Wyświetlenie wyników
    fprintf('FAR: %.4f\n', FAR);
    fprintf('FRR: %.4f\n', FRR);
    fprintf('Współczynnik dyskryminacji D: %.4f\n', D);

    % Macierz pomyłek
    TP = sum(genuine_scores < threshold);
    FN = sum(genuine_scores >= threshold);
    FP = sum(impostor_scores < threshold);
    TN = sum(impostor_scores >= threshold);

    confusion_matrix = [TP, FP; FN, TN];
    disp('Macierz pomyłek:');
    disp(confusion_matrix);

    % Krzywa DET
    [FPR, TPR, ~, AUC] = perfcurve([ones(1, length(genuine_scores)), zeros(1, length(impostor_scores))], ...
                                   [genuine_scores, impostor_scores], 1);
    figure;
    plot(FPR, 1 - TPR);
    xlabel('False Positive Rate (FAR)');
    ylabel('False Negative Rate (FRR)');
    title('Krzywa DET');

    % Rozkłady wyników
    figure;
    histogram(genuine_scores, 'Normalization', 'pdf', 'DisplayName', 'Genuine Scores');
    hold on;
    histogram(impostor_scores, 'Normalization', 'pdf', 'DisplayName', 'Impostor Scores');
    legend;
    title('Rozkłady wyników osób przyjętych i odrzuconych');
end

function distance = hammingDistance(code1, code2)
    if length(code1) ~= length(code2)
        error('Kody muszą być tej samej długości');
    end
    distance = sum(code1 ~= code2) / length(code1); % Normalizowana odległość Hamminga
end
