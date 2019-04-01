%CODE FOR THE SOFTWARE OF 3D MAPPING SYSTEM
>> release(colorDevice);   %release the color device from the Kinect in order to be able to aquire new one
>> release(depthDevice);   %release the depth device from the Kinect in order to be able to aquire new one
>> imaqhwinfo; %get the info about the center
>> colorDevice = imaq.VideoDevice('kinect', 1);    %activate the colorDevice to aquire an image
>> depthDevice = imaq.VideoDevice('kinect', 2);    %activate the depthDevice to aquire an image
>> colorImage10 = step(colorDevice);   %get RGB image to be processed
>> depthImage10 = step(depthDevice);   %get depth image to be processed
>> release(colorDevice);   %release the color device from the Kinect in order to stop aquiring
>> release(depthDevice);   %release the depth device from the Kinect in order to stop aquiring
>> colorImage = colorImage015;  %create object to hold the RGB image
>> depthImage = depthImage015;  %create object to hold the Depth Image
>> surface = surf(depthImage);  %create surface object to represent the Depth Data
>> XData_depthImage = Surface.XData;    %retrieve X-data for the depth surface
>> YData_depthImage = Surface.YData;    %retrieve Y-data for the depth surface
>> ZData_depthImage = Surface.ZData;    %retrieve Z-data for the depth surface
>> colorImage_cropped = imcrop (colorImage, [700 200 499 699]); %cropped the RGB image firstly at initial 500x700px to get rid of the       
                                                                %background noise and stuff
>> depthSurface = surf(depthImage); %save the surface object at more proper name
>> Xmax-depth = max(XData_depthImage);  %retrieve the maximum X-value (it will represent the pixel size of the image)
>> Ymax-depth = max(YData_depthImage);  %retrieve the maximum Y-value (it will represent the pixel size of the image)
>> Central_Value_Depth = ZData_depthImage((Ymax-depth / 2), (Xmax-depth / 2));  %Calculate the central depth, because we consider this to
                                                                                %be the most accurate way to determine the distance from
                                                                                %the Kinect to the ground
>> Color_crop_gray = rgb2gray(colorImage_cropped);  %get the Grayscale version of the cropped RGBpped_RGB
>> Color_crop_BW=im2bw(Color_crop_gray,0.25);   %create Black&White image using the threshold 0.25
>> Color_crop_Imcomplement = imcomplement (Color_crop_BW);  %revert the B&W image, because we will examine the black circle dots 
                                                            %and it is better for the MATLAB to see them as white ones
>> Color_crop_STREL=imclose(Color_crop_Imcomplement, strel('disk',8)); %Use disk STREL filter with area 8px to fill in the blanks 
>> Labels_Color_crop = bwlabel(Color_crop_STREL);   %Create labels for the detected objects
>> AreaObj_Color_crop = []; %create empty vector to hold areas
>> PerimeterObj_Color_crop = []; %create empty vector to hold perimeters
>> FormatFactorObj_Color_crop = []; %create empty vector to hold formfactors
>> Centroid_Color_crop = [];    %create empty vector to hold centroid coordinates


***************AREAS/ASTRONOMY********************


>> Stats_Color_crop = regionprops (Labels_Color_crop, 'Area', 'Centroid', 'Perimeter'); %get the certain data
                                                                             %for every object
>> for k = 1:numel(Stats_Color_crop) %start a loop to go through every object
    AreaObj_Color_crop = [AreaObj_Color_crop Stats_Color_crop(k).Area]; %Assign area of current object k
                                                    %add object area to the register of all objects
    PerimeterObj_Color_crop = [PerimeterObj_Color_crop Stats_Color_crop(k).Perimeter]; %Assign perimeter of current object k
                                                        %add perimeter to the register of all objects
    FormatFactorObj_Color_crop = [FormatFactorObj_Color_crop ((4*pi*(Stats_Color_crop(k).Area))/((Stats_Color_crop(k).Perimeter)^2))];
                                                        %add form factor to the register of all objects
    Centroid_Color_crop = [Centroid_Color_crop Stats_Color_crop(k).Centroid]; %Assign centroid of current object k
                                                    %add object area to the register of all objects
    end
    
                        %subplot the reaults till now to see if the  automatic selection of borders was correct
>> subplot (2,3,1);     %subplot the images from the steps
>> imshow (colorImage015);       %plot image 
>> title ('RGB');       %set name
>> subplot (2,3,2);     %subplot the images from the steps
>> imshow (colorImage_cropped);       %plot image
>> title ('Cropped_RGB');     %set name
>> subplot (2,3,3);     %subplot the images from the steps
>> imshow (Color_crop_gray);       %plot image
>> title ('Grayscale');     %set name
>> subplot (2,3,4);     %subplot the images from the steps
>> imshow (Color_crop_BW);       %plot image
>> title ('Threshold');      %set name
>> subplot (2,3,5);     %subplot the images from the steps
>> imshow (Color_crop_Imcomplement);       %plot image
>> title ('Complement');        %set name
>> subplot (2,3,6);     %subplot the images from the steps
>> imshow (Color_crop_STREL);       %plot image
>> title ('Strell-filter'); %set name

%*************FASTER ALTERNAIVE************

>> Color_crop_Xmin_frame = 700; %it will hold the minimum value of X_coord amongst the 4 border black points
                                %it is set to maximum in order for us to be able to reduce it
>> Color_crop_Xmax_frame = 0;   %it will hold the maximum value of X_coord amongst the 4 border black points
                                %it is set to minimum in order for us to be able to increase it
>> Color_crop_Ymin_frame = 500; %it will hold the minimum value of Y_coord amongst the 4 border black points
                                %it is set to maximum in order for us to be able to reduce it
>> Color_crop_Ymax_frame = 0;   %it will hold the maximum value of Y_coord amongst the 4 border black points
                                %it is set to minimum in order for us to be able to reduce it
>> for k = 1:numel(Stats_Color_crop) %start a loop to go through every object
    if (Stats_Color_crop(k).Area > 360 && Stats_Color_crop(k).Area < 435 && ((4*pi*(Stats_Color_crop(k).Area))/((Stats_Color_crop(k).Perimeter)^2)) > 1)
                                            %the main IF selects only te objects that have proper area and their for factor is near/above
                                            % 1 meaning they are circles
        if(Centroid_Color_crop(2*k) < Color_crop_Xmin_frame)    %compare with the current min X value
            Color_crop_Xmin_frame = Centroid_Color_crop(2*k);
        end
        if(Centroid_Color_crop(2*k) > Color_crop_Xmax_frame)    %compare with the current max X value
            Color_crop_Xmax_frame = Centroid_Color_crop(2*k);
        end
        if(Centroid_Color_crop((2*k)-1) < Color_crop_Ymin_frame)    %compare with the current min Y value
            Color_crop_Ymin_frame = Centroid_Color_crop((2*k)-1);
        end
        if(Centroid_Color_crop((2*k)-1) > Color_crop_Ymax_frame)    %compare with the current mmax Y value
            Color_crop_Ymax_frame = Centroid_Color_crop((2*k)-1);
        end
    end
end
>> final_color_crop = imcrop (colorImage_cropped, [Color_crop_Ymin_frame Color_crop_Xmin_frame (Color_crop_Ymax_frame - Color_crop_Ymin_frame) (Color_crop_Xmax_frame - Color_crop_Xmin_frame)]);   %crop the final Color image based on the borders
                                                                            %calculated by the black circle markers
>> delta_X_col_frame = (Color_crop_Xmax_frame - Color_crop_Xmin_frame);     %find the  whole dimention of the image on the X axis
>> delta_Y_col_frame = (Color_crop_Ymax_frame - Color_crop_Ymin_frame);     %find the whole dimention of the image on the Y axis
>> Center_depth_X = 256;    %locate the X_coord of the center of the depth image
>> Center_depth_Y = 212;    %locate the Y_coord of the center of the depth image
>> X_min_depth_crop = (Center_depth_X - (delta_X_col_frame/2)); %calculate the X values corresponding to the same section in the depth image as it is in the color
%we tried to easily allign the both (color and depth) images but that turned out not being so easy
>> Normalized_flipped_surf = ZData_depthImage;  %the folowing operations will normalized and flip the depth image in order
                                                %for us (and for the robot) to be abla to clearly distinguish the buildings
                                                %first we make a copy of the Zdata for the depth Image
>> for i = 1:Xmax-depth                     %first we go through every pixel both on the X and Y axis
        for j = 1:Ymax-depth
            if Normalized_flipped_surf(j,i) > Central_Value_Depth   %here we eliminate the errors from the Kinect itself
                                                                    %Sometimes values above the max distance occur 
                Normalized_flipped_surf(j,i) = Central_Value_Depth; %we equalize them to the maximum distance on order later their to 
                                                                    %be flat when we flip the image
            end
            if Normalized_flipped_surf(j,i) <= 2000                 %also because the Kinect calculates the distance in mm and our sensor 
                                                                    %is placed over 2 meters above the buldings and just like real-life
                                                                    %satelites buildings are much smaller i our case bellow 30cm
                                                                    %thus every value that is too close to the Kinect is considered
                                                                    %as false and nonrelevant therefore reduced value
                Normalized_flipped_surf(j,i) = Central_Value_Depth; %again we equalize to the max distance to make it flat later
            end
            if Normalized_flipped_surf(j,i) == 0                    %at this point value 0 means that the object just next to the 
                                                                    %Kinect sensor which is impossible, thus these values area 
                                                                    %also normalized
                Normalized_flipped_surf(j,i) = Central_Value_Depth;
            end
            Normalized_flipped_surf(j,i) = Central_Value_Depth - Normalized_flipped_surf(j,i);  %with this equation we flip the image
                                                                        % 180 degrees to be better understandable which is what
        end
    end
>> for i = 1:140                                    %the following for-loops compensate for the both left and right areas 
                                                    %where our camera caught part from the dest at which we have established our sensor
        for j = 1:Ymax-depth
            Normalized_flipped_surf(j,i) = 0;
        end
    end
>> for i = 370:512                                  %this does the  same as the previous function
        for j = 1:Ymax-depth
            Normalized_flipped_surf(j,i) = 0;
        end
    end
>> for i = 1:Xmax-depth                             %with this loop we remove and flatten the road there are minor inequalities
        for j = 1:Ymax-depth
            if Normalized_flipped_surf(j,i) <= 150  %the value is about 30% less than the size of the buildings therefore only
                                                    %the real buildings remain on the #d surface image
                Normalized_flipped_surf(j,i) = 0;
            end
        end
    end
>> W/o_airbags_BW_map = imbinarize(Normalized_flipped_surf);    %we binarze the image making it Black&White in order to be able to 
                                                                %make a graph and to calculate the fastest routes
>> W/o_airbags_Complement_BW_map = imcomplement(BW_map);        %we make the complement image in order to turn the 'roads' white 
                                                                %because the function we use to make the graph conects all the WHITE 
                                                                %pixels to their neighbourhood WHITE pixels
                                                                %Thus all the roads are conected in a big Graph
>> for i = 16:497                                %But before doing so we will add the so called 'air bags' which will compensate 
                                                %for the real size of the robot (car) 
        for j = 16:409
            if Normalized_flipped_surf(j,i) > 10
                for k = (j-15):(j+15)               %the value of 15px is derived from the size of the black markers from the 
                                                    %colored image from where we derived the scale [px/mm]
                    for m = (i-15):(i+15)
                        if Normalized_flipped_surf(k,m) < 10
                            Normalized_flipped_surf(k,m) = 10;  %we do not need a big value for the air bags
                        end
                    end
                end
            end
        end
    end
>> surf (XData_depthImage, YData_depthImage, Normalized_flipped_surf);  %we plot the resulting surface to check the result
>> BW_map = imbinarize(Normalized_flipped_surf);        %we make it a binary (B&W) image 
>> Complement_BW_map = imcomplement(BW_map);            %also complementing it
>> Graph = binaryImageGraph(Complement_BW_map, 4);  %here wee create the graph using the MATLAB toolbox - 'Image Graphs'
>> Cropped_binary_final = imcrop (Complement_BW_map, [150 0 300 350]);  %here we crop the final surface in order to concentrate mostly on our field and disgard the desk and other obstacles that got into our field of view
>> Cropped_Graph = binaryImageGraph(Cropped_binary_final, 4);   %to save calculation time we build graph that is not connected to its diagonal neighbours
%we know the coordinates of our 'cities' on the image
%they are A(150,125), O(150,300), C(60,125), K(50,300)
>> h_for_graph_plot = plot(Cropped_Graph,'Layout','force'); %this is used for the highlighting function

>> G = Cropped_Graph; %create object to hold the final graph
>> [rows,cols] = size(G.Nodes); %get the size of the Nodes data
>> for m = 1:rows   %go through every node
       if (G.Nodes{m,1} == 150 && G.Nodes{m,2} == 125)  %if the coordinated coincide one of the cities we take its pixel Index
            A_px_Index = G.Nodes{m,3};
       end
       if (G.Nodes{m,1} == 150 && G.Nodes{m,2} == 300)
            O_px_Index = G.Nodes{m,3};
       end
       if (G.Nodes{m,1} == 60 && G.Nodes{m,2} == 125)
            C_px_Index = G.Nodes{m,3};
       end
       if (G.Nodes{m,1} == 50 && G.Nodes{m,2} == 300)
            K_px_Index = G.Nodes{m,3};
       end
   end
>> pathA_O = shortestpath(G,A_px_Index,O_px_Index); %use the function shortest path which in this case derive the path using
                                                    %the  Dijkstra`s algrithm for fastest route
>> highlight(h_for_graph_plot,pathA_O,'NodeColor','g','EdgeColor','g'); %we highlight the paths on the final graph
>> Sub = subgraph(G, pathA_O);  %we can create a subgraph to hold the paths for every couple of cities
>>  plotImageGraph(G)       %finaly we plot subgraph above the main one for representation of the desired path
    hold on
    plotImageGraph(Sub)