function evaluate_iris_model(image_folder, codes_folder, threshold)
    % image_folder: folder z obrazami tęczówek
    % codes_folder: folder z kodami tęczówek
    % threshold: próg decyzyjny dla przyznania dostępu

    num_persons = 21; % Liczba osób w bazie danych
    images_per_person = 3; % Liczba zdjęć na osobę

    % Tablice do przechowywania wyników
    genuine_scores = [];
    impostor_scores = [];

    % Ocena dla każdej osoby
    for i = 1:num_persons
        for j = 1:images_per_person
            image_path = fullfile(image_folder, num2str(i), sprintf('image%d.jpg', j));
            new_iris_code = iris_recognition_ofta(image_path);
            
            for k = 1:num_persons
                for l = 1:images_per_person
                    code_path = fullfile(codes_folder, num2str(k), sprintf('code%d.mat', l));
                    load(code_path, 'iris_code');
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
    plot(FPR, TPR);
    xlabel('False Positive Rate (FAR)');
    ylabel('True Positive Rate (1 - FRR)');
    title('Krzywa DET');

    % Rozkłady wyników
    figure;
    histogram(genuine_scores, 'Normalization', 'pdf', 'DisplayName', 'Genuine Scores');
    hold on;
    histogram(impostor_scores, 'Normalization', 'pdf', 'DisplayName', 'Impostor Scores');
    legend;
    title('Rozkłady wyników osób przyjętych i odrzuconych');
end