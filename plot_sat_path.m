function [] = plot_sat_path( lat_array, lon_array )
%% This function reads in a set of lattitude and Longitude arrays
% it plots them onto a project map
% Input Vector format
%   lat_array(number of passes, data for each pass);
%   lon_array(number of passes, data for each pass); 

%% generates a world map with color coded elevation and continets
map = worldmap([-89.90 89.90],[-180 180]);
% framem('off');
gridm('off');
mlabel('off');
plabel('off');
tightmap
load topo
load geoid
geoshow(topo, geoidrefvec, 'DisplayType', 'texturemap');
land=shaperead('landareas','UseGeoCoords',true);
geoshow([land.Lat], [land.Lon], 'Color','k');    

%% generatates test orbits
% domain = -180:1:180;
% x=ones(10,max(size(domain)));
% y=ones(10,max(size(domain)));
% step = (97/60)*24;
% for i=1:size(x,1)
%     for j=1:max(size(x))
%         x(i,j) = domain(1,j);
%           if((i==3)&&(j<100)&&(j>75))
%             y(i,j) = NaN;
% %             x(i,j) = NaN;
%           else
%             y(i,j) = 40.5*sin( (x(i,j)-(i-1)*step)*pi/180 );
%           end
%     end
% end
% lat_array = x;
% lon_array = y;

%% Places satellite lat and lon on the map
lat_size = size(lat_array);
lat = zeros(1,lat_size(2));
lon = zeros(1,lat_size(2));
for k=1:lat_size(1)
    for i=1:lat_size(2)
        if( ~isnan(lat_array(k,i)) )
            lat(1,i) = lat_array(k,i);
        else
            lon(1,i) = NaN;
            lat(1,i) = NaN;
        end

        if( ~isnan(lon_array(k,i)) ) 
            lon(1,i) = lon_array(k,i);
        else
            lon(1,i) = NaN;
            lat(1,i) = NaN;
        end
    end
    geoshow(lon,lat,'Color','w','DisplayType','line');
end


end