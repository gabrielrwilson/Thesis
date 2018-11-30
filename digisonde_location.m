function [number, name, latitude, longitude, date_first_data, date_last_data] = digisonde_location_giro()
    %% Open the Digisonde locations document
    filename = 'Digisonde_names_locations_giro';
    [num, txt] = xlsread(filename);

    %% Get the number and locations of the digisonde

    [n,~]=size(num);
    number = num(2:n,1);
%     URSI = txt(2:n,2);
    name = txt(2:n,3);
    latitude = num(2:n,3);
    longitude = num(2:n,4);
    date_first_data_raw = txt(2:n,6);
    date_last_data_raw = txt(2:n,7);

    %% Reformat the dates
    date_first_data = char(n-1,11);
    date_last_data = char(n-1,11);
    for i=1:(n-1)
        first = char(date_first_data_raw(i,1));
        last = char(date_last_data_raw(i,1));

        k_first = 1;
        k_last = 1;
        first_stop = 0;
        last_stop = 0;
        
        [~,first_size] = size(first);
        for j=1:first_size
            if( (j~=first_size)&&(j>1) )
                if(first(j+1)=='(')
                    first_stop = 1;
                elseif(first(j-1)==',')
                    first_stop = 0;
                end     
            end
            if(first_stop == 0)
                date_first_data(i,k_first) = first(1,j);
                k_first = k_first+1;
            end
        end
        
        [~,last_size] = size(last);
        for j=1:last_size
            if( (j~=first_size)&&(j>1) )
                if(last(j+1)=='(')
                    last_stop = 1;
                elseif(last(j-1)==',')
                    last_stop = 0;
                end
            end
            if(last_stop == 0)
                date_last_data(i,k_last) = last(j);
                k_last = k_last+1;
            end
        end        
        datenum(date_first_data, 'mmm dd yyyy')
    end

end

