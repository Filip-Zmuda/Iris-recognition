    function distance = hammingDistance(code1, code2)
        if length(code1) ~= length(code2)
            error('Kody muszą być tej samej długości');
        end
        distance = sum(code1 ~= code2) / length(code1); % Normalizowana odległość Hamminga
    end