function [Df] = calcula_Df(w, f, d, Vprop, barrido)

    N = 7;
    Df = zeros(length(f),length(barrido));

    for freq=1:length(f)
        for i=1:length(barrido)
            Dsum = 0;
            for n=1:N
                Dsum = Dsum + exp(2 * pi * 1j * f(freq) * n * d * cos(barrido(i)) / Vprop) * conj(w(freq, n));

            Df(freq, i) = Dsum;
            end
        end
    end

    Df = abs(Df);

end