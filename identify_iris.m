% Funkcja główna do porównania kodu tęczówki nowego zdjęcia z bazą danych
function [matched_person, minHammingDistance, access_granted] = identify_iris(image_path, codes_path, threshold)
    % Oblicz kod tęczówki dla nowego zdjęcia
    new_iris_code = iris_recognition_ofta(image_path);

    % Inicjalizacja minimalnej odległości Hamminga i identyfikatora osoby
    minHammingDistance = inf;
    matched_person = -1;

    % Przeszukiwanie bazy danych kodów tęczówek
    for i = 1:21
        person_folder = fullfile(codes_path, num2str(i));
        code_files = dir(fullfile(person_folder, '*.mat'));
        
        for j = 1:length(code_files)
            % Załaduj kod tęczówki z bazy danych
            load(fullfile(person_folder, code_files(j).name), 'iris_code');
            
            % Oblicz odległość Hamminga
            hammingDist = hammingDistance(new_iris_code, iris_code);
            
            % Aktualizacja minimalnej odległości Hamminga i identyfikatora osoby
            if hammingDist < minHammingDistance
                minHammingDistance = hammingDist;
                matched_person = i;
            end
        end
    end
    
    % Decyzja o przyznaniu dostępu
    access_granted = minHammingDistance < threshold;

    % Wyświetl wynik
    if access_granted
        fprintf('Dostęp przyznany. Najbardziej podobne oko należy do osoby nr %d.\n', matched_person);
        fprintf('Odległość = %f.\n', minHammingDistance);
    else
        fprintf('Odmowa dostępu. Najbardziej podobne oko należy do osoby nr %d.\n', matched_person);
        fprintf('Minimalna odległość = %f.\n', minHammingDistance);
    end
end
