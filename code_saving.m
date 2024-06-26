clear all;
close all;

% Ścieżki folderów
base_path = 'C:\Users\filip\Documents\studia dell\6 semestr\Biometria\OFTA\OFTA_reshape\pictures';
save_path = 'C:\Users\filip\Documents\studia dell\6 semestr\Biometria\OFTA\OFTA_reshape\codes';

% Sprawdzenie, czy folder 'codes' istnieje, jeśli nie, tworzenie go
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Przeglądanie folderów ponumerowanych od 1 do 21
for i = 1:21
    % Ścieżka do bieżącego folderu
    current_folder = fullfile(base_path, num2str(i));
    
    % Ścieżka do zapisu kodów tęczówek
    save_folder = fullfile(save_path, num2str(i));
    
    % Sprawdzenie, czy folder zapisu istnieje, jeśli nie, tworzenie go
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    
    % Pobieranie listy plików BMP w bieżącym folderze
    image_files = dir(fullfile(current_folder, '*.bmp'));
    
    % Przetwarzanie każdego zdjęcia w folderze
    for j = 1:length(image_files)
        % Ścieżka do bieżącego zdjęcia
        image_path = fullfile(current_folder, image_files(j).name);
        
        % Generowanie kodu tęczówki
        iris_code = iris_recognition_ofta(image_path);
        
        % Ścieżka do zapisu kodu tęczówki
        save_file = fullfile(save_folder, [image_files(j).name(1:end-4) '.mat']);
        
        % Zapisanie kodu tęczówki do pliku .mat
        save(save_file, 'iris_code');
    end
end
