function tswrite(fileName, freq, Y)
% tswrite(fileName, freq, Y)
%
% Writes multiport network parameters in Touchstone format
% Now only handles y-parameters with up to 4 ports
%

comment = 'no comment supplied';

f = fopen(fileName, 'wt');

fprintf(f, '%s\n', ['! ' comment]);
fprintf(f, '%s\n', '# Hz Y RI R 50');

N = size(Y, 1);
nf = size(Y, 3);
for i = 1:nf
    fprintf(f, '%.11e', freq(i));
    y = Y(:,:,i);
    if N <= 2
        % 1 or 2 ports - column-major, entire matrix per line
        for nn = 1:(N*N)
            fprintf(f, ' %.11e %.11e', real(y(nn)), imag(y(nn)));
        end
        fprintf(f, '\n');
    else
        % 3+ ports - matrix row per line
        offset = '';
        for m = 1:N
            fprintf(f, '%s', offset);
            for n = 1:N
                fprintf(f, ' %.11e %.11e', real(y(m,n)), imag(y(m,n)));
            end
            fprintf(f, '\n');
            offset = '               '; % for second and subsequent rows
        end
    end
end

fclose(f);
