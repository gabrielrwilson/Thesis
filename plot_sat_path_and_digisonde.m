function [overlap_lat,overlap_lon,overlap_alt,overlap_time,overlap_name,overlap_distance,overlap_ursi] = plot_sat_path_and_digisonde(lat_array,lon_array,alt_array,time_array,map_name)
%% This function reads in a set of lattitude and Longitude arrays
% it plots them onto a project map
% Input Vector format
%   lat_array(number of passes, data for each pass);
%   lon_array(number of passes, data for each pass); 

%% generates a world map with color coded elevation and continets
h1 = worldmap([-89.90 89.90],[-180 180]);
title_start = strfind(map_name,'STPSat3_');
map_title = map_name(title_start:title_start+16);
map_title = strrep(map_title,'STPSat3_','');
map_title = strcat(map_title(1:2),'_',map_title(3:5),',_',map_title(6:9));
map_title = strrep(map_title,'_', ' ');


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
title({map_title, 'Digisonde Overflights'});
%% generatates test orbits
% domain = -180:1:180;
% x=ones(10,max(size(domain)));
% y=ones(10,max(size(domain)));
% step = (97/60)*24;
% for i=1:size(x,1)
%     for j=1:max(size(x))
%         x(i,j) = domain(1,j);
%           if((i==3)&&(j<100)&&(j>75))
% %             y(i,j) = NaN;
% %             x(i,j) = NaN;
%           else
%             y(i,j) = 40.5*sin( (x(i,j)-(i-1)*step)*pi/180 );
%           end
%     end
% end
% lon_array = x;
% lat_array = y;

%% Places satellite lat and lon on the map
lat_size = max(size(lat_array));
k_number=min(size(lat_array));
sat_lat = zeros(1,lat_size);
sat_lon = zeros(1,lat_size);
sat_alt = zeros(1,lat_size);
sat_time = zeros(1,lat_size);
for k=1:k_number
    for i=1:lat_size
        if( ~isnan(lat_array(k+i-1)) )
            sat_lat(1,i) = lat_array(k+i-1);
            sat_alt(1,i) = alt_array(k+i-1);
            sat_time(1,i)=time_array(k+i-1);
        else
            sat_lon(1,i) = NaN;
            sat_lat(1,i) = NaN;
        end

        if( ~isnan(lon_array(k+i-1)) ) 
            sat_lon(1,i) = lon_array(k+i-1);
        else
            sat_lon(1,i) = NaN;
            sat_lat(1,i) = NaN;
        end
    end
    h1 = geoshow(sat_lat,sat_lon,'Color','w','DisplayType','line','DisplayName','STPSat-3 Path');
end    
%% Place Digisonde locations on map
[ ds_ursi, ds_name, ds_lat, ds_lon,~,~] = digisonde_location_giro();
ds_lat = ds_lat';
ds_lon = ds_lon';
ds_name= ds_name';

lat_size = size(ds_lat);
ds1_lat = zeros(1,lat_size(2));
ds1_lon = zeros(1,lat_size(2));
ds1_name = cell(1,lat_size(2));
% ds1_ursi = cell(1,lat_size(2));
for k=1:lat_size(1)
    for i=1:lat_size(2)
        drop_loc = 0;

        if( ~isnan(ds_lat(k,i)) && (drop_loc ==0) )
            ds1_lat(1,i) = ds_lat(k,i);
        else
            ds1_name(1,i) = 'none';
%             ds1_usri(1,i) = 'none';
            ds1_lon(1,i) = NaN;
            ds1_lat(1,i) = NaN;
            drop_loc = 1;
        end

        if( ~isnan(ds_lon(k,i)) && (drop_loc ==0) )
            if(ds_lon(k,i)>180)
                ds1_lon(1,i) = ds_lon(k,i)-360;
            else
                ds1_lon(1,i) = ds_lon(k,i);
            end
            ds1_name(1,i) = ds_name(k,i);
%             ds1_usri(1,i) = ds_usri(k,i);
        else
            ds1_name(1,i) = 'none';
%             ds1_usri(1,i) = 'none';
            ds1_lon(1,i) = NaN;
            ds1_lat(1,i) = NaN;
        end           
    end
    h2 = geoshow(ds1_lat,ds1_lon,'Color','w','DisplayType','point','DisplayName','Digisonde');
end
% textm(ds1_lat,ds1_lon,ds1_name);

%% Determine where the satellite path passes over a Digisonde
sat_lat_size = max(size(sat_lat));
k_number = min(size(sat_lat));
ds_lat_size = size(ds_lat);

calibration_distance = 140;  %300 km radius of
% calibration_deg = km2deg(calibration_distance,6878.1);

% calibration_distance = 1000;  %300 km radius of
overlap_lat1 = zeros(1,sat_lat_size);
overlap_lon1 = zeros(1,sat_lat_size);
overlap_alt1 = zeros(1,sat_lat_size);
overlap_time1 = zeros(1,sat_lat_size);
arc_distance = zeros(1,sat_lat_size);
overlap_name1= cell(1,sat_lat_size);
overlap_ds = nan(2,sat_lat_size);
overlap_ursi1= cell(1,sat_lat_size);

m=1;
for k=1:k_number
for j=1:ds_lat_size(2)
%     for i=1:max(sat_lat_size)
        arc_length = distance('gc',sat_lat(k,:),sat_lon(k,:),ds_lat(1,j),ds_lon(1,j),(6378.1+sat_alt(k,i)));        
        access = find(arc_length < calibration_distance);
        if(min(size(access))>0)
            for l=1:max(size(access))         
                overlap_lat1(1,m) = sat_lat(k,access(l));
                overlap_lon1(1,m) = sat_lon(k,access(l));
                overlap_alt1(1,m) = sat_alt(k,access(l));
                overlap_time1(1,m) = sat_time(k,access(l));
                overlap_name1(1,m)= ds_name(1,j);
                overlap_ursi1(1,m)= ds_ursi(j,1);
                overlap_ds(1,m)= ds_lat(1,j);
                overlap_ds(2,m)= ds_lon(1,j);
                arc_distance(1,m) = arc_length(access(l));
                m=m+1;  
%                 break;
            end
        end
%     end
end
end

%% Place Digisonde passes on map
if(m>1)
    access_number = 0;
    for i=1:(max(size(ds_name))-1)
        access_name = ds_name(1,i);
        if(strcmp(access_name,ds_name(1,i+1)) )
           access_number = access_number+1;
           
           
        end
    end    
    overlap_lat = overlap_lat1( 1,1:(m-1) );
    overlap_lon = overlap_lon1( 1,1:(m-1) );
    overlap_alt = overlap_alt1( 1,1:(m-1) );
    overlap_time = overlap_time1( 1,1:(m-1) );
    overlap_name = overlap_name1( 1,1:(m-1) );
    overlap_ursi = overlap_ursi1(1,1:(m-1) );

    overlap_distance = arc_distance( 1,1:(m-1) );
    h3 = geoshow(overlap_ds(1,:),overlap_ds(2,:),'Marker','o','MarkerFaceColor','w','Color','k','DisplayType','point');
    legend([h1 h2 h3],'STPSat-3 Path','Digisonde','Possible Calibration','Location','SouthEast');
else
    overlap_lat = 0;
    overlap_lon = 0;
    overlap_alt = 0;
    overlap_time = 0;
    overlap_name = 0;
    
    legend([h1 h2],'STPSat-3 Path','Digisonde','Location','SouthEast');
end

% textm(overlap_lat,overlap_lon+1,overlap_name,'FontSize',8,'Rotation',15);    
% plotname = strrep(sourceFileName,'L3.nc','');
% plotname=strcat(plotname,'map_',num2str(i));
    poutfile=['print(''-djpeg'',''',map_name,''')'];
    eval(poutfile)
end